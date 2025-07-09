# Week 1: Intro to Unix and the Command Line

Welcome to Week 1! 
This week, you'll get comfortable using the Unix command line, which is an essential skill for working with genomic data on any server or HPC.

---
## 🧠 Learning Objectives

By the end of this exercise, you should be able to:

- Navigate the file system using `cd`, `ls`, and `pwd`
- Create, move, copy, and delete files and directories
- Use wildcards and pipes
- Read and manipulate files using `head`, `tail`, `less`, `grep`, `cut`, and `sort`
- Understand and run a basic bash script

---
## 🧪 Exercise 1: Accessing LEAP2 - TXST's high-performance computer (HPC)
Connecting to an HPC system is typically done with a program known as “SSH” (Secure SHell) which is accessed through a Terminal. Linux and Mac users will find a Terminal program already installed on their computers. Windows users will need to install a Terminal emulator. I suggest using MobaXterm. If you do not already have it installed, please install the free version now: https://mobaxterm.mobatek.net/download.html. 

The SSH program needs the web address of the server we want to access. You should have received an email from Shane Flaherty with your login information and the LEAP2 web address.
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


## 🧪 Exercise 2: Navigating the Filesystem

Let's check where you are by <ins>p</ins>rinting your <ins>w</ins>orking <ins>d</ins>irectory (i.e. where you currently are in the system).
```bash
pwd
```

Now <ins>l</ins>i<ins>s</ins>t the contents of your current directory using the `ls` command.
```bash
ls
```

Sometimes there are *hidden* files, typically configuration files, which begin with a dot (`.`). E.g., your bash profile is configured by the file ~/.bash_profile. To see these files use the **'ls'** *command* with the *argument* -a, to see <ins>a</ins>ll files, including hidden ones.
```bash
ls -a
```

We can also use `ls` to see the sizes of the files, in bytes, in our directories with the argument -l.
```bash
ls -l
```

### I prefer to add -lh, the "h" prints the sizes in a <ins>h</ins>uman-readable format.
```bash
ls -lh
```

Ok, lets move into other directories, or <ins>c</ins>hange <ins>d</ins>irectories using the `cd` command. Move to the XX directory and use `pwd` and `ls` to see where you are and what is there.
```bash
cd xx/
pwd
ls
```

To move back one directory, simply use `..`:
```bash
cd ..
```

You can also use `cd` to "jump" to other directories quickly, like below:
> [!TIP]
> If we are trying to specify a file or path we can begin typing its name and then press the <ins>tab</ins> key to complete it (try it out below). If there is only one possible way to finish what we’ve started typing, it will complete it entirely for us. If there is more than one possible way to finish what we’ve started typing, it will complete as far as it can, and then hitting tab twice quickly will show all the possible options. If tab-complete does not do either of those things, then we are either confused about where we are or what is where, or we've maybe spelled the name wrong.

```bash
cd xx/xxx/xxxx/
pwd
```

When combined with "..", you can move back multiple directories as well.
```bash
cd ../xx/xxx/xxxx/
pwd
```

Lastly, to return all the way back to your root directory, simple just use `cd`.
```bash
cd 
pwd
```

## 🧪 Exercise 3: Creating, copying, moving, and removing files + directories
> [!WARNING]
> Using commands that do things like create, copy, and move files at the command line will <ins>**overwrite**</ins> files if they have the same name. And using commands that delete things will do so <ins>**permanently**</ins>. Use caution using these commands.

Start with <ins>m</ins>a<ins>k</ins>ing a new <ins>dir</ins>ectory, using the `mkdir` command, then move into it.
```bash
mkdir practice_files
cd practice_files
```

Now create files in this directory and use `ls` to check that the files are there.
```bash
echo "Sample 1" > sample1.txt
echo "Sample 2" > sample2.txt
echo "Sample 3" > sample3.txt
ls
```

Lets <ins>c</ins>o<ins>p</ins>y the second file and rename it. The first file called is the file you want to copy, the second is what you want to name the copy.
```bash
cp sample2.txt copy_sample2.txt
ls
```

To make a copy and put it somewhere else, like in new subdirectory, we can change the second positional argument using a relative path (“relative” because it starts from where we currently are). Make a new directory called "copies" then copy "Sample 3" there.
```bash
mkdir copies
cp sample3.txt copies/copy_sample3.txt
cd copies/
ls
```

Ok, now we want to move the copy_sample2.txt file here in the copies/ directory. To do that we will use the <ins>m</ins>o<ins>v</ins>e command `mv`), but remember the file is in the previous directory (hint: ".."). 
```bash
mv ../copy_sample2.txt .
ls
```

Here, I am also introducing the use of a single `.`, which is telling `mv` that I want to move the file to where I currently am. 

The `mv` command can also be used to rename files like so:
```bash
mv copy_sample2.txt renamed_copy_sample2.txt
ls
```

Lets get rid of the renamed file. To <ins>r</ins>e<ins>m</ins>ove files, use the `rm` command. 
```bash
rm renamed_copy_sample2.txt
ls
```

You can also remove entire directories. Let's remove the copies/ directory using -r. 
```bash
cd ..
rm -r copies/
ls
```

## 🧪 Exercise 4: Editing/creating file contents

It is often very useful to be able to generate new plain-text files quickly at the command line, or make some changes to an existing one. One way to do this is using a text editor that operates at the command line. Here we’re going to look at one program that does this called `nano`. Let's test it with a file that already exists.
```bash
nano sample1.txt
```

This will open up an interface and allow you to add text. Add two sample names, A and B.
```bash
sampleA
sampleB
```

To save the file and exit, we need to use some of the keyboard shortcuts listed on the bottom. Type "ctrl" + "x". It will ask if you want to save, type "y" and then press "enter". Alternatively, if you wanted to change the file name, you can before pressing enter. 

To get a "sneak-peak" at the file, we can either use the `head` command to show the top of the file contents. There is also `tail` , which prints the last 10 lines of a file by default:
```bash
head sample1.txt
```

If we wanted to count the number of lines, words, or characters a file has, we can use the `wc` or (<ins>w</ins>ord <ins>c</ins>ount) command. 
```bash
wc sample1.txt
```

To *only* get the number of lines in the file, use the argument -l.
```bash
wc -l sample1.txt
```

I also use nano to create files, you just simply type `nano`, followed by the file name you want and its extension. E.g.:
```bash
nano mynewfile.txt
```

>[!NOTE]
> The second part of a file name is called the filename extension, and indicates what type of data the file holds. Here are some common examples:
>* .txt is a plain text file.
>* .csv is a text file with tabular data where each column is separated by a comma.
>* .tsv is like a CSV but values are separated by a tab.
>* .log is a text file containing messages produced by a software while it runs.
>* .pdf indicates a PDF document.
>* .png is a PNG image.
>* .sh indicates a shell script file.
>  
>This is just a convention, albeit an important one. Files contain bytes: it’s up to us and our programs to interpret those bytes according to the rules for plain text files, PDF documents, configuration files, images, and so on.
>Naming a PNG image of a whale as whale.mp3 doesn’t somehow magically turn it into a recording of whalesong, though it might cause the operating system to try to open it with a music player when someone double-clicks it.


## 🧪 Exercise 5: Redirectors and Wildcards

Go back to XX directory. 

Redirectors change where the output of a command is going. Lets look at an example using a "pipe" (`|`).
Remember our commands `ls` and `wc -l`? If we “pipe” (`|`) the `ls` command into the `wc -l` command, instead of printing the output from `ls` to the screen as usual, it will go into `wc -l` which will print out how many items there are in your current working directory:
```bash
ls | wc -l
```

Another important character is the greater than sign, `>`. This tells the command line to redirect the output to a file, rather than just printing it to the screen as we’ve seen so far. For an example of this we will write the output of ls to a new file called “directory_contents.txt”:
```bash
ls > directory_contents.txt
head directory_contents.txt
```
**It’s important to remember that the `>` redirector will overwrite the file we are pointing to if it already exists.**
If we want to append output to a file, rather than overwrite it, we can use two signs, `>>` instead.

Now, wildcards are special characters that enable us to specify multiple items very easily. Let’s say we only want to look for files, in our current working directory, that end with the extension “.txt”. The * wildcard can help us with that.
```bash
ls *.txt
```
What this is saying is that no matter what comes before, if it ends with “.txt” we want it.

Test: move all .txt files into a new directrory called new_dir
