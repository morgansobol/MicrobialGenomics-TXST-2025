curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
press Enter
type yes, press Enter
Press Enter to confirm location 
Do you wish to update your shell profile to automatically initialize conda?
type Yes, press Enter

Close out of terminal window and reopen, or source ~/.zshrc to restart
zaj33 >> to (base) zaj33 --> signifies the base conda environment

conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 

