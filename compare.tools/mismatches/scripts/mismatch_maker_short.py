import random
#A script that takes in reads and puts mismatches in them
#Example: Input - 5 reads with no mismatches
# Output - 5 reads with 1,2,3,4,5 mismatches 
with open("short.reads.mismatches.fastq") as f:
    lines = f.readlines()

basepairs = "ATCG"

#Loop through each read
for i in range(1,len(lines),4):
    numLoops = (i-1)/4+1
    #Create numLoops mismatches in read
    alreadyReplaced = [False]*len(lines[i])
    for j in range(int(numLoops)):
        #TODO: figure out if position has already been replaced
        while True:
            position = random.randint(0,len(lines[i])-2)
            if (alreadyReplaced[position] == False):
                alreadyReplaced[position] = True
                break
        while True:
            replace = random.choice(basepairs)
            if (replace != lines[i][position]):
                break
        lines[i] = lines[i][:position] + replace + lines[i][position+1:]
        readNameLength=len(lines[i-1])
    #Change read names
    lines[i-1]= lines[i-1][:readNameLength-1] + str(numLoops) + lines[i-1][readNameLength-1]

with open("short.reads.mismatches.fastq",'w') as f:
    f.writelines(lines)

