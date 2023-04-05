for folder in ["sparse", "dense","super-dense"]:
    for i in range(1,26):
        filename = "./{}/{}task.task".format(folder,i)  
        print(filename)
        with open(filename, 'r') as f:
            lines = f.readlines()
            num_lines = len(lines)
        with open(filename, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(str(num_lines) + '\n' + content)
