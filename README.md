# How to run:

requires ffmpeg for producing the videos

## Cartoon of pinch (in drift)

```
cd Cartoon_drift
python 000_buildup.py
cd Cartoon
python 001_pinch.py
cd ../..
./copy_frames.sh
```

Slide-by-slide frames can then be found in folder named frame_by_frame_cartoon

## Pinch video (in drift)

Instead of Drift, the same applies for Dipole and Quadrupole.

```
cd Drift
python 000_buildup.py
cd Pinch
python 001_pinch.py
python 002_frame_pinch.py
chmod +x 003_make_movie.sh
./003_make_movie.sh
```

Pinch movie will be then available as Drift/Pinch/pinch.mp4

