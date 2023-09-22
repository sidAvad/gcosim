#Python script takes in a bpfile and modifies it by introducing gcos accordign to a user specified list of 
#simulation parameters 
#TODO:Add documentation everywhere
#TODO:Add assertions for number of chromosomes etc.
import numpy as np
import pandas as pd 
import sys
from bisect import bisect
import math,random
from datetime import datetime
#test edit remode copy
#test of github deskop workflow

def sample_breakpoints(mapDF_c,sex,options):

    physical_positions=mapDF_c['pos'].tolist()

    if sex == 'male':
        genetic_positions=mapDF_c['male_cM'].tolist()
    elif sex == 'female':
        genetic_positions=mapDF_c['female_cM'].tolist()
    else:
        sys.exit('invalid sex specification')   

    centimorgan_range = [genetic_positions[0] , genetic_positions[-1]]

    #Add independent gcos
    gco_midpoints_genetic=sample_exponentialRecombinations(centimorgan_range,options['indep_gco_rate'])
    gco_midpoints_physical = np.interp(gco_midpoints_genetic,genetic_positions,physical_positions)
    gco_breakpoints_physical=[x for y in [[(math.floor(x - np.rint(np.random.normal(loc=options['length']/2, scale=options['length']/50))),False),(math.floor(x + np.rint(np.random.normal(loc=options['length']/2, scale=options['length']/50))),True)] for x in gco_midpoints_physical ] for x in y ]
    #print('indep gco bps : {}'.format(gco_breakpoints_physical))
    
    recomb_breakpoints_genetic=sample_exponentialRecombinations(centimorgan_range)
    recomb_breakpoints_physical = [(math.floor(x),False) for x in np.interp(recomb_breakpoints_genetic,genetic_positions,physical_positions)] #(note the use of math.floor to round down physical positions from floats to ints )

    all_breakpoints_physical=recomb_breakpoints_physical+gco_breakpoints_physical

    #Add complex gcos
    complexgco_breakpoints=[]   
    for (bp,gco_status) in recomb_breakpoints_physical: #TODO: Bug : We should be looping only over the recombination breakpoints.
        dice_roll=random.choice(list(range(options['beta'])))
        if dice_roll == 1:
            complexgco_loc=bp + np.rint(np.random.normal(loc=0, scale=options['dist']))
            complexgco_breakpoints.append((complexgco_loc-np.rint(np.random.normal(loc=options['length']/2, scale=options['length']/50)),False)) #TODO:use truncated normal (or maybe a different distribution)
            complexgco_breakpoints.append((complexgco_loc+np.rint(np.random.normal(loc=options['length']/2, scale=options['length']/50)),True))
        else:
            continue
    
    #print('complexgco bps : {}'.format(complexgco_breakpoints))
    breakpoints_final=all_breakpoints_physical+complexgco_breakpoints
    breakpoints_final.sort(key=lambda y:y[0])

    return(breakpoints_final)

def sample_exponentialRecombinations(centimorgan_range,rate=1):
    '''
    Sample successive points along a genetic_map given the map range and the exponential rate parameter (rate in morgan units)
    '''
    
    scale = 100/rate #Adjustment for centimorgan scale of input
    start = centimorgan_range[0]
    end = centimorgan_range[1]

    breakpoints = np.array([]) #Initialize empty array to store breakpoints
    sampler_current_position = start

    while True: #while loop to track whether sampler has reached end of sample
        next_position = sampler_current_position + np.random.exponential(scale)
        if next_position > end:
            break
        breakpoints = np.append(breakpoints,next_position)
        sampler_current_position = next_position

    return(breakpoints)

def transmit_segments(breakpoints,inherited_segments):
    '''
    generate transmitted segment give a pair of inherited segments and a set of breakpoints
    '''
    if len(breakpoints)==0:
        writesegs=random.choice(inherited_segments) #Pick haplotype at random when there are no recombinations 
    else:
        writesegs=[]
        inherited_segmentbounds=[]
        for h in [0,1]:
            inherited_segmentbounds.append([d['start'] for d in inherited_segments[h]] + [inherited_segments[h][-1]['end']])
        j=1 #index for inherited_segments
        #copying_hap=0 #start copying path from haplotype 0 for each chromosome; TODO: random(0,1) 
        copying_hap=random.randint(0,1) #start copying path randomly from haplotype 0 or haplotype 1 for each chromosome
        writeseg=dict(inherited_segments[copying_hap][j-1])
        for bpno,(r,gco_endpoint) in enumerate(breakpoints):
            #Add first set of segments (added segments are always behind current recombination)
            while (writeseg['start'] < r):
                if writeseg['end'] < r:
                    #print('added segments : {};{}'.format(writeseg,j))
                    writesegs.append(writeseg)
                    j=j+1
                    writeseg=dict(inherited_segments[copying_hap][j-1])
                else:
                    writeseg['end']=r
                    if (gco_endpoint==True): #Check if gco_endpoint is true
                        writeseg['gco']=True
                    else: #If not , segment is a non gco segment
                        writeseg['gco']=False
                    #print('added segments and switch : {};{}'.format(writeseg,j))
                    writesegs.append(writeseg)
                    prevhap=writeseg['hap'] # record the previous haplotype broken by bp as recipient haplotype
                    copying_hap=not copying_hap  # Switch haplotype only if we reach the recombination
                    j=bisect(inherited_segmentbounds[copying_hap],r) #Find segment recombination falls in on switched haplotype
                    writeseg=dict(inherited_segments[copying_hap][j-1])
                    writeseg['start']=r
                    writeseg['hapopp']=prevhap # record the recipient haplotype of current haplotype
                    break

        #Add segments after final recombination 
        while (j <= len(inherited_segments[copying_hap])):
            writeseg=dict(inherited_segments[copying_hap][j-1])
            writeseg['start']=r
            writeseg['gco']=False #Final segment cannot be gco because it is bounded by chromosome end
            writesegs.append(writeseg)
            j+=1

    return(writesegs)

def merge_adjacent(inherited_segments_hap):
    merged_segments_hap=[]
    prev_d={'start':None,'stop':None,'hap':None, 'hapopp':None} #HAPOPP
    for i,d in enumerate(inherited_segments_hap):
        if d['hap']==prev_d['hap']:
            new_d=d
            new_d['start']=prev_d['start']
            prev_d=new_d
        else:
            merged_segments_hap.append(prev_d)
            prev_d=d

    merged_segments_hap.append(prev_d)

    merged_segments_hap_final=[d for d in merged_segments_hap if not d['hap']==None]
    return(merged_segments_hap_final)

class _Node():

    def __init__(self,idno,sex,mapDF):
        
        self.idno=idno
        self.sex=sex
        self.mapDF=mapDF
        self.gco_params={'dist':5000,'length':1000,'beta':5,'indep_gco_rate':10}
        self.bpstr_recomb=[[],[]]
        self.df_gco=None

        self.breakpoints=[[]]*22
        self.transmitted_segments=[[None]]*22
        self.chromosome_endpoints=[(mapDF[mapDF["#chr"]==c].iloc[0].pos, mapDF[mapDF["#chr"]==c].iloc[-1].pos) for c in range(1,23)]

        assert (sex == "male") | (sex == "female")
        assert self.idno > 0

    def generate_breakpoints(self):
        
        for c in range(1,23):
            mapDF_c=self.mapDF[self.mapDF["#chr"]==c]

            bp_in_map_range=False
            while not bp_in_map_range: #Resample breakpoints until they fall inside the range of mapfile
                self.breakpoints[c-1]=sample_breakpoints(mapDF_c,self.sex,self.gco_params)
                if all([x[0] < mapDF_c.iloc[-1].pos for x in self.breakpoints[c-1]]):
                    bp_in_map_range=True

    def transmit_segments(self):
        '''
        generate transmitted segment give a pair of haplotype ids and a set of breakpoints
        '''
        for c in range(1,23):
            self.transmitted_segments[c-1]=transmit_segments(self.breakpoints[c-1],self.inherited_segments[c-1])

    def merge_segments(self):
        '''
        Pull out recombination segments and merge them
        '''

        for h in [0,1]:
            for c in range(1,23):
                self.merged_recomb_segments[c-1][h]=merge_adjacent([d for d in self.inherited_segments[c-1][h] if d['gco']==False])

    def generate_bplines_recomb(self):

        headers=['g{}_i{}_s{}_h0'.format(self.generation,self.idno,self.simno),'g{}_i{}_s{}_h1'.format(self.generation,self.idno,self.simno)]
    
        for h in [0,1]:
            #Generate recomb lines only if segments separated by recombinations have been merged
            assert not self.merged_recomb_segments==None
            bpstr_list=[]
            for c in range(1,23):
                l=self.merged_recomb_segments[c-1][h]
                str_c = '{}|'.format(c)+str(int(l[0]['start']))+' '+''.join(['{}:{} '.format(d['hap'],int(d['end'])) for d in l if d['gco']==False])[:-1] #TODO:Convert to integer in transmit_segments rather than here
                bpstr_list.append(str_c + ' ')
                self.bpstr_recomb[h]=headers[h] + ' ' + ''.join(bpstr_list)[:-1]
    
    def generate_bplines_gco(self):
        headers=['g{}_i{}_s{}_h0'.format(self.generation,self.idno,self.simno),'g{}_i{}_s{}_h1'.format(self.generation,self.idno,self.simno)]
        gco_dict_list = []
        for h in [0,1]:
            #gco_table=pd.DataFrame()
            
            for c in range(1,23):
                l=self.inherited_segments[c-1][h]
                #str_c = '{}|'.format(c)+str(int(l[0]['start']))+' '+''.join(['{}/{}:{} '.format(d['hap'],int(d['start']),int(d['end'])) for d in l if d['gco']==True])[:-1] #TODO:Convert to integer in transmit_segments rather than here
                #gco_table.append(str_c + ' ')
                #self.bpstr_gco[h]=headers[h] + ' ' + ''.join(gco_table)[:-1]
                for d in l:
                    if d['gco']:
                        d['chr'] = c
                        gco_dict_list.append(d)
        self.df_gco=pd.DataFrame(gco_dict_list)

class ChildNode(_Node):
    '''
    '''
    def __init__(self,idno,sex,mapDF,ParentNodes):
        assert len(ParentNodes)==2
        assert not ParentNodes[0].sex==ParentNodes[1].sex
        assert not ParentNodes[0].idno==ParentNodes[1].idno
        assert ParentNodes[0].simno==ParentNodes[1].simno
        assert ParentNodes[0].generation==ParentNodes[1].generation
        
        super().__init__(idno,sex,mapDF)
        self.generation = ParentNodes[0].generation + 1
        self.simno=ParentNodes[0].simno
        
        tmp_segments=[[],[]]
        #Transmit segments here and update attributes for inherited_haplotype
        for i,parent in enumerate(ParentNodes):
            parent.generate_breakpoints()
            parent.transmit_segments()
            tmp_segments[i]=parent.transmitted_segments

        self.inherited_segments=[list(c) for c in zip(*tmp_segments)]
        self.merged_recomb_segments=[[None,None] for i in range(22)]


class FounderNode(_Node):
    '''
    Datastructure that holds founder individuals 
    '''
    def __init__(self,idno,sex,simno,mapDF):
        super().__init__(idno,sex,mapDF)
        self.simno=simno
        self.inherited_segments=[\
                #[[{'start':self.chromosome_endpoints[c-1][0],'end':self.chromosome_endpoints[c-1][1],'hap':str(2*self.idno-1),'gco':False}],\
                #[{'start':self.chromosome_endpoints[c-1][0],'end':self.chromosome_endpoints[c-1][1],'hap':str(2*self.idno),'gco':False}]]\
                [[{'start':self.chromosome_endpoints[c-1][0],'end':self.chromosome_endpoints[c-1][1],'hap':str(2*self.idno-1),'hapopp':str(2*self.idno),'gco':False}],\
                [{'start':self.chromosome_endpoints[c-1][0],'end':self.chromosome_endpoints[c-1][1],'hap':str(2*self.idno),'hapopp':str(2*self.idno-1),'gco':False}]]\
                for c in range(1,23)]
        self.generation=0


if __name__ == "__main__":

    def simulate_admixed_file(T,nsim,outfile,mapDF):

        t=T
        tree=[[] for _ in range(T+1)]
        for i in range(2**T):
            if i%2:
                tree[T-t].append(FounderNode(i+1,'male',nsim,mapDF))
            else:
                tree[T-t].append(FounderNode(i+1,'female',nsim,mapDF))
        
        t=t-1

        while t >= 0:
            print(t)
            for i in range(2**t):
                if i%2:
                    print("creating male : i == {}",format(i))
                    tree[T-t].append(ChildNode(i+1,'male',mapDF,[tree[T-t-1][2*i],tree[T-t-1][2*i+1]]))
                else:
                    print("creating female : i == {}",format(i))
                    tree[T-t].append(ChildNode(i+1,'female',mapDF,[tree[T-t-1][2*i],tree[T-t-1][2*i+1]]))
            t=t-1
        
        tree[T][0].merge_segments() 
        tree[T][0].generate_bplines_recomb()
        tree[T][0].generate_bplines_gco()


        with open(outfile+'_{}.recomb.bp'.format(nsim),'w') as f:
            f.write('\n'.join(tree[T][0].bpstr_recomb))

        #with open(outfile+'_{}.gco.bp'.format(nsim),'w') as f:
        #    f.write('\n'.join(tree[T][0].df_gco))
        tree[T][0].df_gco.to_csv(outfile+'_{}.gco.bp'.format(nsim), sep='\t', index=False)

        #print(tree[T][0].inherited_segments)

        return()

    T=int(sys.argv[1])
    nsim=int(sys.argv[2])
    mapDF=pd.read_table(sys.argv[3])
    
    seed = datetime.now().timestamp()
    random.seed(seed)
    print("Seed used:", seed)
    
    simulate_admixed_file(T,nsim,'gcos',mapDF)
