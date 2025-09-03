## For Mac users
```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
```
```bash
bash Miniconda3-latest-MacOSX-arm64.sh
```
Type yes, press Enter.

Press Enter again to confirm location.

Do you wish to update your shell profile to automatically initialize conda?
Type Yes, press Enter. 

Close out of terminal window and reopen, or run below to restart:
```bash
source ~/.zshrc
```

zaj33 >> to (base) zaj33 --> signifies the base conda environment

If you get an error about accepting terms of service, type the following commands. 
```bash
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
```
```bash
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 
```
