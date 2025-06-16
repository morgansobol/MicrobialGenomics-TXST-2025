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

### 🔹 Let's check where you are by <ins>P</ins>rinting your <ins>W</ins>orking <ins>D</ins>irectory (i.e. where you currently are in the system)
```bash
pwd
```

### 🔹 Now <ins>L</ins>i<ins>S</ins>t the contents of your current directory
```bash
ls
```

### 🔹 Sometimes there are *hidden* files, typically configuration files, which begin with ".". For example, your bash profile is configured with ~/.bash_profile. To see these files use ls, but with -a, to see <ins>A</ins>ll files, incliding hidden ones
```bash
ls -a
```

### 🔹 We can also use ls to see the sizes of the files in our directories by appending -l
```bash
ls -l
```

### But I prefer to add -lh, the "h" prints the sizes in a <ins>H</ins>uman-readable format
```bash
ls -lh
```

### 🔹 Ok, lets move into other directories, or <ins>C</ins>hange <ins>D</ins>irectories using the "cd" command. Move to the XX directory and use "pwd" and "ls" to see where you are and what is there
```bash
cd xx/
pwd
ls
```

### 🔹 To move back one directory, simply do:
```bash
cd ..
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

### 🔹 <ins>C<i/ns>o<ins>P</ins>y the second file and rename it. The first file called is the file you want to copy, the second is what you want to name the copy.
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

## 🧪 Exercise 3: Using wildcards


