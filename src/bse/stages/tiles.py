"""Stage `tiles`: one 224-px tile at 0.5 µm/px centred on every spot of every section with a
full-resolution H&E image (DESIGN §6), and the frozen encoders' CLS features, cached per section.

The spot centre in full-resolution pixels is the Space Ranger pixel position scaled by the ratio of
the full-resolution image width to the Space Ranger input image width (recovered from the hires scale
factor); the tile side in full-resolution pixels is 112 µm / (µm per pixel). Features: Phikon (ViT-B/16,
iBOT, through the Hugging Face hub) and an ImageNet ViT-B/16 (timm). Both read the same tiles.
"""
import json
import os
import sys
import time

THREADS = int(os.environ.get("BSE_THREADS", "4"))
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, str(THREADS))           # one thread pool, before numpy and torch load

import numpy as np
import pandas as pd
from PIL import Image

from bse import config as C
from bse import provenance as P
from bse import sections as S

Image.MAX_IMAGE_PIXELS = None
OUT = C.RESULTS / "tiles"
FEAT = C.DATA / "features"
BATCH = 32
SPOT_UM = 55.0                            # a Visium spot is 55 µm across; Space Ranger records its pixel diameter
TILE_UM = C.TILE_PX * C.TILE_MPP          # 112 µm


def load_image(path) -> np.ndarray:
    path = str(path)
    if path.lower().endswith((".tif", ".tiff")):
        import tifffile
        with tifffile.TiffFile(path) as tf:
            page = tf.pages[0]
            arr = page.asarray()
    else:
        arr = np.asarray(Image.open(path).convert("RGB"))
    if arr.ndim == 2:
        arr = np.stack([arr] * 3, axis=-1)
    if arr.shape[-1] > 3:
        arr = arr[..., :3]
    return arr


def spot_centres(section: S.Section, img_shape) -> tuple[np.ndarray, float, float]:
    """Full-resolution pixel centres (x, y), the scale from Space Ranger pixels to full resolution, and
    the full-resolution µm per pixel (from the 55 µm spot diameter Space Ranger recorded)."""
    sf = S.scalefactors(section)
    hires_long = 2000                                 # Space Ranger's hires image is 2000 px on its long side
    input_long = hires_long / sf["tissue_hires_scalef"]
    scale = max(img_shape[:2]) / input_long
    mpp_full = (SPOT_UM / sf["spot_diameter_fullres"]) / scale
    x = section.adata.obs.pxl_col.values * scale
    y = section.adata.obs.pxl_row.values * scale
    return np.column_stack([x, y]), scale, mpp_full


def extract_tiles(img: np.ndarray, centres: np.ndarray, side_px: int) -> np.ndarray:
    h, w = img.shape[:2]
    half = side_px // 2
    tiles = np.zeros((len(centres), C.TILE_PX, C.TILE_PX, 3), dtype=np.uint8)
    for i, (x, y) in enumerate(centres):
        x0, y0 = int(round(x)) - half, int(round(y)) - half
        x1, y1 = x0 + side_px, y0 + side_px
        sx0, sy0, sx1, sy1 = max(x0, 0), max(y0, 0), min(x1, w), min(y1, h)
        patch = np.full((side_px, side_px, 3), 255, dtype=np.uint8)       # white where the tile leaves the image
        if sx1 > sx0 and sy1 > sy0:
            patch[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = img[sy0:sy1, sx0:sx1]
        tiles[i] = np.asarray(Image.fromarray(patch).resize((C.TILE_PX, C.TILE_PX), Image.BILINEAR))
    return tiles


class Encoders:
    def __init__(self):
        import timm
        import torch
        from transformers import AutoModel
        torch.set_num_threads(THREADS)
        self.torch = torch
        self.phikon = AutoModel.from_pretrained("owkin/phikon").eval()
        self.imagenet = timm.create_model(C.ENCODERS["imagenet"], pretrained=True, num_classes=0).eval()
        self.mean = torch.tensor([0.485, 0.456, 0.406]).view(1, 3, 1, 1)
        self.std = torch.tensor([0.229, 0.224, 0.225]).view(1, 3, 1, 1)

    def features(self, tiles: np.ndarray) -> dict[str, np.ndarray]:
        torch = self.torch
        out = {"phikon": [], "imagenet": []}
        with torch.no_grad():
            for i in range(0, len(tiles), BATCH):
                x = torch.from_numpy(tiles[i:i + BATCH]).permute(0, 3, 1, 2).float() / 255.0
                x = (x - self.mean) / self.std
                out["phikon"].append(self.phikon(pixel_values=x).last_hidden_state[:, 0].numpy())
                out["imagenet"].append(self.imagenet(x).numpy())
        return {k: np.concatenate(v) for k, v in out.items()}


def run():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    FEAT.mkdir(parents=True, exist_ok=True)
    secs = [s for s in S.load_all() if s.image is not None]
    print(f"  {len(secs)} sections with a full-resolution image")
    enc = None
    rows = []
    for s in secs:
        done = FEAT / f"{s.id}.npz"
        if done.exists():
            z = np.load(done, allow_pickle=True)
            rows.append({"section": s.id, "dataset": s.dataset, "subtype": s.subtype, "spots": len(z["barcodes"]),
                         "scale": float(z["scale"]), "mpp_full_res": float(z["mpp"]), "side_px": int(z["side_px"]),
                         "seconds": float(z["seconds"])})
            print(f"  skip    {s.id} (features cached)"); continue
        t1 = time.time()
        img = load_image(s.image)
        print(f"    {s.id}: image {img.shape[1]}x{img.shape[0]} decoded in {time.time() - t1:,.0f}s", flush=True)
        centres, scale, mpp = spot_centres(s, img.shape)
        side_px = int(round(TILE_UM / mpp))
        tiles = extract_tiles(img, centres, side_px)
        del img
        print(f"    {s.id}: {len(tiles)} tiles extracted at {time.time() - t1:,.0f}s; encoding", flush=True)
        if enc is None:
            enc = Encoders()
        feats = enc.features(tiles)
        # a thumbnail strip of 8 tiles per section, for the page and for checking alignment by eye
        strip = np.concatenate(list(tiles[np.linspace(0, len(tiles) - 1, 8).astype(int)]), axis=1)
        Image.fromarray(strip).save(OUT / f"tiles_{s.id}.png")
        np.savez_compressed(done, barcodes=s.adata.obs_names.values.astype(str), phikon=feats["phikon"],
                            imagenet=feats["imagenet"], centres=centres, scale=scale, side_px=side_px,
                            mpp=mpp, seconds=time.time() - t1)
        rows.append({"section": s.id, "dataset": s.dataset, "subtype": s.subtype, "spots": len(tiles),
                     "scale": scale, "mpp_full_res": mpp, "side_px": side_px, "seconds": round(time.time() - t1, 1)})
        print(f"  {s.id:<16} {len(tiles):>5} tiles, side {side_px}px ({mpp:.3f} µm/px), scale {scale:.3f}, "
              f"{time.time() - t1:,.0f}s", flush=True)
    inv = pd.DataFrame(rows)
    inv.to_csv(OUT / "tiles.csv", index=False)
    (OUT / "encoders.json").write_text(json.dumps({"phikon": "owkin/phikon (ViT-B/16, iBOT) via the Hugging Face hub",
                                                   "imagenet": C.ENCODERS["imagenet"], "tile_px": C.TILE_PX,
                                                   "tile_mpp": C.TILE_MPP, "threads": THREADS}, indent=1))
    P.log_run("tiles", {"tile_px": C.TILE_PX, "tile_mpp": C.TILE_MPP, "threads": THREADS},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*"))], time.time() - t0)
    print(f"done in {(time.time() - t0) / 60:.0f} min")


if __name__ == "__main__":
    sys.exit(run())
