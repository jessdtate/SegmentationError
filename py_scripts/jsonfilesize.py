import sys, json

usage  = "extracts field value from json string \n   Usage: python jsonfilesize.py json_string field_name"

def main():

    n = len(sys.argv)
        
    if n==1:
        print(usage)
        return
    elif n<3:
        raise ValueError("not enough inputs\n "+usage)
    elif n>3:
        raise ValueError("unexpected inputs\n "+usage)
        
    print(json.loads(sys.argv[1])[sys.argv[2]])
    
    return



if __name__ == "__main__":
    main()
