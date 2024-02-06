
import sys,os,subprocess
from subprocess import PIPE
def main(argv):

    if "-h" in argv or "?" in argv:
        print("This script interprets a qstat call and prints with better formatting/control. The width of the jobname field can be set with the -w argument and jobs from a specific user can be requested using the -u argument")
        quit()
        
    # Save defaults
    user=""
    name_width=100
    while len(argv) > 0:
        if argv[0] == '-w':
            argv.pop(0)
            if len(argv) == 0:
                print("ERROR: expected an integer to follow the -w flag. Exiting...")
                quit()
            else:
                name_width=int(argv.pop(0))
        elif argv[0] == "-u":
            argv.pop(0)
            if len(argv) == 0:
                print("ERROR: expected a user name to follow the -u flag. Exiting...")
                quit()
            else:
                user=argv.pop(0)
        else:
            print("ERROR: could not interpret the argument {}. Exiting...".format(argv[0]))
            quit()


    # redirect a qstat call into output
    p = subprocess.Popen("qstat -f".split(),stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate(b"input data that is passed to subprocess' stdin")
    output = str(output)

    # Initialize job information dictionary
    dict = {}
    for count_i,i in enumerate(output.split(r'\n')):
        fields = i.split()
        if len(fields) == 0: continue
        if "Job Id" in i:
            current_key = i.split()[2]
            dict[current_key] = { "State":"NA" , "Name":"NA", "Walltime":"NA", "Queue":"NA", "User":"NA"}
            continue
        if "Job_Name" == fields[0]:
            dict[current_key]["Name"] = fields[2]
        if "job_state" == fields[0]:
            dict[current_key]["State"] = fields[2]
        if "queue" == fields[0]:
            dict[current_key]["Queue"] = fields[2]
        if "Resource_List.walltime" == fields[0]:
            dict[current_key]["Walltime"] = fields[2]        
        if "Job_Owner" == fields[0]:
            dict[current_key]["User"] = fields[2].split("@")[0]

    template_string="{:<40s} {:<"+str(name_width)+"s} {:<20s} {:20s} {:<20s}"
#    print "{:<40s} {:<100s} {:<20s} {:20s} {:<20s}".format("ID","Name","User","Walltime","Status","Queue")
    print(template_string.format("ID","Name","User","Walltime","Status","Queue"))
    for i in list(dict.keys()):
        if user == "" or dict[i]["User"] == user:
            print(template_string.format(i,dict[i]["Name"],dict[i]["User"],dict[i]["Walltime"],dict[i]["State"],dict[i]["Queue"]))

if __name__ == "__main__":
   main(sys.argv[1:])
