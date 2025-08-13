## 🧪 Exercise 1: Accessing LEAP2 - TXST's high-performance computer (HPC)
Connecting to an HPC system is typically done with a program known as “SSH” (Secure SHell) which is accessed through a Terminal.

Linux and Mac users will find a Terminal program already installed on their computers. Windows users will need to install a Terminal emulator. I suggest using MobaXterm. If you do not already have it installed, please install the free version now: https://mobaxterm.mobatek.net/download.html. 

Now, the _Terminal_ is a text input and output environment where we can type commands and see the output. In other words, it is the "window" in which you enter the actual commands and those commands are interpreted and run by a _Shell_. 

So the _Shell_ is the program inside the terminal that actually processes commands and returns the output. In most Linux and Mac operating systems, it uses a _Bash_ shell, which is essentially its own programming language and what we will use below. 

In summary, think of it this way: Terminal is the TV and Shell is the program running the TV. 

The SSH program allows us to connect to a remote computer/server over a network and execute commands and transfer files to this remote server. So, to access TXST's remote HPC, SSH requires the web address of the server. You should have received an email from Shane Flaherty with your login information and the LEAP2 web address.

Open your terminal to test connecting to LEAP2. 

```bash
ssh [insert net ID]@leap2.txstate.edu
```
If this is your first time connecting to LEAP2, you will probably see a similar message as the one below. 
```bash
The authenticity of host 'leap2.txstate.edu' can't be established.
RSA key fingerprint is 2a:b6:f6:8d:9d:c2:f8:2b:8c:c5:03:06:a0:f8:59:12.
Are you sure you want to continue connecting (yes/no)?
```
This is your computer warning you that you’re about to connect to another computer, type "yes" and press ENTER to proceed. This will add the HPC to your "known hosts", and you shouldn’t see the message again the future.

You should now be prompted to input your password. Your first password can be found in the email from Shane. Type it in carefully because no characters will appear on the screen for you to see what you type. Then press ENTER.

If you entered your password appropriately, congratulations, you’re now connected to the HPC!

Now, if it is your first time accessing LEAP2, we will change your password as per Shane's request.
Type `passwd' and press ENTER to start 
