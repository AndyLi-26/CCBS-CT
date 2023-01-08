import os
if __name__=="__main__":
    a=os.listdir()
    for i in a:
        if i[-3:]!="csv":
            continue
        print(i)
        with open(i,'r') as f:
            content=[]
            for l in f:
                l=l.strip().strip(',')
                if l=='':
                    break
                l=l[l.index("tasks")+3:]
                l=l[l.index("tasks")+5:]
                l=l.replace(".xml",'')
                content.append(l)
            print(content)
            
        with open(i,'w') as f:
           print(*content,sep='\n',file=f) 
    
    
