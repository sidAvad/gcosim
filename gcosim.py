#Python script takes in a bpfile and modifies it by introducing gcos accordign to a user specified list of 
#simulation parameters 
import numpy as np
import pandas as pd 
import sys
from bisect import bisect
import math,random
from datetime import datetime

def sample_breakpoints(mapDF_c,sex,options):

    physical_positions=mapDF_c['pos'].tolist()

    if sex == 'male':
        genetic_positions=mapDF_c['male_cM'].tolist()
    elif sex == 'female':
        genetic_positions=mapDF_c['female_cM'].tolist()
    else:
        sys.exit('invalid sex specification')   

    centimorgan_range = [genetic_positions[0] , genetic_positions[-1]]
    # options['length']: fixed gco length    
    # TODO: Fixed the bug of gco breakpoint samping. Oct 24, 2023
    
    # 1. Sampling recombination breakpoints
    recomb_breakpoints_genetic = sample_exponentialRecombinations(centimorgan_range)
    recomb_breakpoints_physical = [(math.floor(x),False,None) for x in np.interp(recomb_breakpoints_genetic,genetic_positions,physical_positions)] #(note the use of math.floor to round down physical positions from floats to ints )

    # 2. Sampling simple gco starting breakpoints
    # uniformed_genetic_positions
    slope = (genetic_positions[-1] - genetic_positions[0]) / (physical_positions[-1] - physical_positions[0])
    intercept = genetic_positions[0] - slope * physical_positions[0]
    uniformed_genetic_positions = slope * np.array(physical_positions) + intercept
    mixed_genetic_positions = options['alpha'] * uniformed_genetic_positions + (1 - options['alpha']) * np.array(genetic_positions)
    
    simplegco_startpoints_genetic = sample_exponentialRecombinations(centimorgan_range, options['simple_gco_rate'])
    simplegco_startpoints_physical = [(math.floor(x), False, 'simple-gco') for x in np.interp(simplegco_startpoints_genetic, mixed_genetic_positions, physical_positions)] # TODO: checking if this data structure makes sence
    
    # 3. Sampling complex gco starting breakpoints
    complexgco_startpoints = []   
    for (bp,gco_status, gco_type) in recomb_breakpoints_physical: 
        dice_roll=random.choice(list(range(options['beta'])))
        if dice_roll == 1:
            complexgco_loc = bp + np.rint(np.random.normal(loc=0, scale=options['dist']))
            complexgco_startpoints.append((math.floor(complexgco_loc - (options['length'] / 2)), False, 'complex-gco'))
        else:
            continue
    
    # 4. Checking for simulated gcos conflicting with other gcos
    gco_startpoints_all = simplegco_startpoints_physical + complexgco_startpoints
    gco_startpoints_all.sort(key = lambda y:y[0])
    is_no_coflict = np.ones(len(gco_startpoints_all))
    # loop for invalid gcos before the starting physical position
    i = 0
    while i < len(gco_startpoints_all) and gco_startpoints_all[i][0] < physical_positions[0]:
        is_no_coflict[i] = 0
        i += 1
    # loop for invalid gcos beyonds the ending physical position
    i = len(gco_startpoints_all) - 1
    while i >= 0 and gco_startpoints_all[i][0] + options['length'] > physical_positions[-1]:
        is_no_coflict[i] = 0
        i -= 1
    # main loop for (1) gcos including a recombination bp; (2) gcos overlapping with previous gcos
    curr_recomb = 0
    for i in range(len(gco_startpoints_all)):
        if is_no_coflict[i]:
            # (1) gcos including a recombination bp
            if len(recomb_breakpoints_physical) > 0 and curr_recomb < len(recomb_breakpoints_physical):
                while curr_recomb < len(recomb_breakpoints_physical) and gco_startpoints_all[i][0] > recomb_breakpoints_physical[curr_recomb][0]:
                    curr_recomb += 1
                if curr_recomb < len(recomb_breakpoints_physical): 
                    if gco_startpoints_all[i][0] <= recomb_breakpoints_physical[curr_recomb][0] and gco_startpoints_all[i][0] + options['length'] >= recomb_breakpoints_physical[curr_recomb][0]:
                        is_no_coflict[i] = 0
                        continue
            # (2) next gcos overlapping with current gco
            j = i + 1
            while j < len(gco_startpoints_all) and gco_startpoints_all[i][0] + options['length'] >= gco_startpoints_all[j][0]:
                is_no_coflict[j] = 0
                j += 1
    gco_startpoints_all_final = [gco_startpoints_all[i] for i in range(len(gco_startpoints_all)) if is_no_coflict[i]]
    gco_endpoints_all_final = [(gco_startpoints_all_final[i][0] + options['length'], True, gco_startpoints_all_final[i][2]) for i in range(len(gco_startpoints_all_final))]    
    
    breakpoints_final = recomb_breakpoints_physical + gco_startpoints_all_final + gco_endpoints_all_final
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
        for bpno,(r,gco_endpoint, gco_type) in enumerate(breakpoints):
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
                        writeseg['gco_type']=gco_type # When the new segment is gco, record its gco type
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

        #Add segments after final recombination # Bugs fixed here
        writeseg=dict(inherited_segments[copying_hap][j-1])
        writeseg['start']=r # update the start point as r
        while (j <= len(inherited_segments[copying_hap])):
            writesegs.append(writeseg)
            if (j == len(inherited_segments[copying_hap])):
                break
            else:
                j=j+1
                writeseg=dict(inherited_segments[copying_hap][j-1])

    return(writesegs)

def merge_adjacent(inherited_segments_hap):
    merged_segments_hap=[]
    prev_d={'start':None,'stop':None,'hap':None, 'hapopp':None}
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
        self.gco_params={'dist':5000,'length':300,'beta':5,'simple_gco_rate':12,'alpha':0.1} # 'alpha': the weight of uniformed genetic map, where 1 - alpha is the weight of genetic map for recombination
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

        headers=['g{}_i{}_h0'.format(self.generation,self.idno),'g{}_i{}_h1'.format(self.generation,self.idno)]
    
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
        headers=['g{}_i{}_h0'.format(self.generation,self.idno),'g{}_i{}_h1'.format(self.generation,self.idno)]
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
                        d['hap_id'] = h
                        gco_dict_list.append(d)
        self.df_gco=pd.DataFrame(gco_dict_list)

class ChildNode(_Node):
    '''
    '''
    def __init__(self,idno,sex,mapDF,ParentNodes):
        assert len(ParentNodes)==2
        assert not ParentNodes[0].sex==ParentNodes[1].sex
        assert not ParentNodes[0].idno==ParentNodes[1].idno
        assert ParentNodes[0].generation==ParentNodes[1].generation
        
        super().__init__(idno,sex,mapDF)
        self.generation = ParentNodes[0].generation + 1
        
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
    def __init__(self,idno,sex,mapDF):
        super().__init__(idno,sex,mapDF)
        self.inherited_segments=[\
                [[{'start':self.chromosome_endpoints[c-1][0],'end':self.chromosome_endpoints[c-1][1],'hap':str(2*self.idno-1),'hapopp':str(2*self.idno),'gco':False,'gco_type':None}],\
                [{'start':self.chromosome_endpoints[c-1][0],'end':self.chromosome_endpoints[c-1][1],'hap':str(2*self.idno),'hapopp':str(2*self.idno-1),'gco':False,'gco_type':None}]]\
                for c in range(1,23)]
        self.generation=0


if __name__ == "__main__":

    def simulate_admixed_file(T,ninds,outprefix,mapDF):

        t=T
        tree=[[] for _ in range(T+1)]
       
        #Create the top layer of 100 founders
        for ind in range(ninds):
            if ind%2:
                tree[T-t].append(FounderNode(ind+1,'male',mapDF))
            else:
                tree[T-t].append(FounderNode(ind+1,'female',mapDF))
        
        t=t-1
        
        
        while t >= 0:
            print(t)
            for ind in range(ninds):
                #select random male and random female founder
                male_founder=random.choice([x for x in tree[T-t-1] if x.sex=='male'])
                female_founder=random.choice([x for x in tree[T-t-1] if x.sex=='female'])
                if ind%2:
                    print("creating male : ind == {}",format(ind))
                    tree[T-t].append(ChildNode(ind+1,'male',mapDF,[male_founder,female_founder]))
                else:
                    print("creating female : ind == {}",format(ind))
                    tree[T-t].append(ChildNode(ind+1,'female',mapDF,[male_founder,female_founder]))
            t=t-1
        
        tree[T][0].merge_segments() 
        tree[T][0].generate_bplines_recomb()
        tree[T][0].generate_bplines_gco()


        with open(outprefix+'.recomb.bp','w') as f:
            for ind in range(ninds):
                f.write('\n'.join(tree[T][ind].bpstr_recomb))
        
        for ind in range(ninds):
            df_gco_output=tree[T][ind].df_gco
            del df_gco_output['gco'] # remove the redundant column 'gco' from output (as they must be gcos)
            tree[T][ind].df_gco.to_csv(outprefix+'.gco.bp', sep='\t', index=False)

        #print(tree[T][0].inherited_segments)

        return()

    T=int(sys.argv[1])
    ninds=int(sys.argv[2])
    outprefix=sys.argv[3]
    mapDF=pd.read_table(sys.argv[4])
    
    seed = datetime.now().timestamp()
    random.seed(seed)
    print("Seed used:", seed)
    
    simulate_admixed_file(T,ninds,outprefix,mapDF)
