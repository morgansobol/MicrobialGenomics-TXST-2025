## For Mac users
```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
```
```bash
bash Miniconda3-latest-MacOSX-arm64.sh
```
Press enter to scroll through the ToS. Then type yes when prompted, press Enter.

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
## For Windows using MobaXTerm

Search for Windows Features and open
Check Windows Subsystem for Linux like below. NOTE: you will need to reboot your computer!
<img width="750" height="665" alt="image" src="https://github.com/user-attachments/assets/d9702436-b41b-47db-9afd-65d87637b199" />

Now search for Windows Powershell and open it. Install Windows Subsystem for Linux (wsl)
```shell
wsl --install
```
Open MobaXTerm, got to Settings, and select the Terminal tab. 
In the Local shell settings box, from the Dropdown select WSL default.
Now your MobaXTerm terminal will always start with WSL. 

Now install Miniconda. 
```shell
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
```bash
bash Miniconda3-latest-MacOSX-arm64.sh
```

Press enter to scroll through the ToS. Then type yes when prompted, press Enter.

Press Enter again to confirm location.

Do you wish to update your shell profile to automatically initialize conda?
Type Yes, press Enter. 

Close out of terminal window and reopen. Conda should now be installed.

zaj33 >> to **(base)** [name] --> signifies the base conda environment






