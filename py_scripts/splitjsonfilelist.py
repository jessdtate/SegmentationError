import sys, json

"""
splitjsonfilelist.py
splits up the output from

rclone lsjson drive:folder

into a list of filenames and sizes.  these can then be parsed with:

list | cut -f1
list | cut -f2

"""

usage  = "python splitjsonfilelist.py json_string"

field_list=["Size","Name"]


def main():

    n = len(sys.argv)
        
    if n==1:
        print(usage)
        return
    elif n>2:
        raise ValueError("unexpected inputs\n "+usage)
    
    contents = json.loads(sys.argv[1])
        
    print(contents)
    print(type(contents))
    print(contents[0])
    print(type(contents[0]))
    print(len(contents))
    
    for c in contents:
#        line=""
#        for f in field_list:
#            line+=c[f]+" "
#        line+="\n"
#    print(line)
        print(str(c["Size"])+" "+c["Name"])
    
    
    return



if __name__ == "__main__":
    main()
