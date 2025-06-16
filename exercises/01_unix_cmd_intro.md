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

## 🧪 Exercise 1: Navigating the Filesystem

### 🔹 Let's check where you are by <ins>P</ins>rinting your <ins>W</ins>orking <ins>D</ins>irectory (i.e. where you currently are in the system).
```bash
pwd
```

### 🔹 Now <ins>L</ins>i<ins>S</ins>t the contents of your current directory.
```bash
ls
```

### 🔹 Sometimes there are *hidden* files, typically configuration files, which begin with dot. E.g., your bash profile is configured by the file ~/.bash_profile. To see these files use ls -a, to see <ins>A</ins>ll files, incliding hidden ones.
```bash
ls -a
```

### 🔹 We can also use ls to see the sizes of the files, in bytes, in our directories by appending -l after ls.
```bash
ls -l
```

### I prefer to add -lh, the "h" prints the sizes in a <ins>H</ins>uman-readable format.
```bash
ls -lh
```

### 🔹 Ok, lets move into other directories, or <ins>C</ins>hange <ins>D</ins>irectories using the "cd" command. Move to the XX directory and use "pwd" and "ls" to see where you are and what is there.
```bash
cd xx/
pwd
ls
```

### 🔹 To move back one directory, simply do:
```bash
cd ..
```

### 🔹 You can also use cd to "jump" to other directories quickly, like below:
> [!TIP]
> If we are trying to specify a file or path we can begin typing its name and then press the <ins>tab</ins> key to complete it (try it out below). If there is only one possible way to finish what we’ve started typing, it will complete it entirely for us. If there is more than one possible way to finish what we’ve started typing, it will complete as far as it can, and then hitting tab twice quickly will show all the possible options. If tab-complete does not do either of those things, then we are either confused about where we are or what is where, or we've maybe spelled the name wrong.

```bash
cd xx/xxx/xxxx/
pwd
```
### 🔹 When combined with "..", you can move back multiple directories as well.
```bash
cd ../xx/xxx/xxxx/
pwd
```
### 🔹 Lastly, to return all the way back to your root directory, simple just use cd.
```bash
cd 
pwd
```

## 🧪 Exercise 2: Creating, copying, moving, and removing files + directories
> [!WARNING]
> Using commands that do things like create, copy, and move files at the command line will <ins>**overwrite**</ins> files if they have the same name. And using commands that delete things will do so <ins>**permanently**</ins>. Use caution using these commands.

### 🔹 <ins>M</ins>a<ins>k</ins>e a new <ins>Dir</ins>ectory and move into it
```bash
mkdir practice_files
cd practice_files
```

### 🔹 Now create files in this directory and use "ls" to check that the files are there
```bash
echo "Sample 1" > sample1.txt
echo "Sample 2" > sample2.txt
echo "Sample 3" > sample3.txt
ls
```

### 🔹 <ins>C</ins>o<ins>P</ins>y the second file and rename it. The first file called is the file you want to copy, the second is what you want to name the copy.
```bash
cp sample2.txt copy_sample2.txt
ls
```

### 🔹 To make a copy and put it somewhere else, like in new subdirectory, we could change the second positional argument using a relative path (“relative” because it starts from where we currently are). Make a new directory called "copies" then copy "Sample 3" there.
```bash
mkdir copies
cp sample3.txt copies/copy_sample3.txt
cd copies/
ls
```
### 🔹 Ok, now we want to move the copy_sample2.txt file here in the copies/ directory. To do that we will use the <ins>m</ins>o<ins>v</ins>e command (mv), but remember the file is in the previous directory (hint: ".."). 
```bash
mv ../copy_sample2.txt .
ls
```
### Here, I am also introducing the use of a single dot, which is telling the mv command that I want to move the file to where I currently am. 

## 🧪 Exercise 3: Using wildcards


