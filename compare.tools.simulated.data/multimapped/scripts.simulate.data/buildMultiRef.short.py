import random
#A script that creates a new reference with the multimapped reads 
with open("shortReads.fq",'r') as f:
    lines = f.readlines()

basepairs = "ATCG"
referenceSequence=""
#testSequence=""
#Loop through each read
for i in range(1,len(lines),4):
    numLoops = (i-1)/4+1
    readLength = len(lines[i])
    #Add each read to referenceSequence depending on read number
            #ie Read 4 will be in the reference sequence 4 times
    for j in range(int(numLoops)):
        #Create random filler sequence (50 bps long)
        randomFiller=""
        for k in range(1,50):
            randomFiller+= random.choice(basepairs)
        referenceSequence+= lines[i][:readLength-1] + randomFiller
#        testSequence+= "A"
with open("multimapped.reference.short.fa",'w+') as f:
    f.writelines(">Reference\n")
    f.writelines(referenceSequence)
    
