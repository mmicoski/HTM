# -*- coding: cp1252 -*-
''' Implementation of the algorithms in HTM_CorticalLearningAlgorithms.pdf

        HIERARCHICAL TEMPORAL MEMORY
        including
        HTM Cortical Learning Algorithms

        VERSION 0.2.1, SEPTEMBER 12, 2011
        ©Numenta, Inc. 2011

    Exctracts of the text are used in comments along the code.        
'''

import random
import math
import pickle
import datetime
import sys
##import profile

############# GENERAL PARAMETERS


### spatial pooler

''' If the permanence value for a synapse is greater than this value,
    it is said to be connected
'''
connectedPerm = 0.3

''' A minimum number of inputs that must be active for a column to be considered during
    the inhibition step.
    IMPORTANT: Needs to be consistent with the expected minimum number of active inputs 
'''
minOverlap = 6

''' A parameter controlling the number of columns that will be winners after the inhibition step '''
desiredLocalActivity = 20

''' Amount permanence values of synapses are incremented during learning. '''
permanenceInc = 0.01


''' Amount permanence values of synapses are decremented during learning. '''
permanenceDec = 0.01


### temporal pooler

''' Number of cells in each column '''
cellsPerColumn = 4

''' Activation threshold for a segment. If the number of active connected synapses in a segment
    is greater than activationThreshold, the segment is said to be active.
'''
activationThreshold = 6


''' The maximum number of synapses added to a segment during learning. '''
newSynapseCount = 4

############# OTHER PARAMETERS

#size of the input vector (x,y)
sinp=(10,10)

#size of the column array (x,y)
scol=(40,50)

#range around connectedPerm for initialization
dperm = 0.1

#number of potential synapses per column
npsyn = 40

#size of region around natural center 
radius = 10

# stored history of overlap and active
historyLen = 1000

# minDutyCycle / maxDutyCycle: hystoryLen should be greater than 1/minDutyFactor
minDutyFactor = 0.01

# maximum boost
maxBoost = 4.0

################################
train = False
logFile = False
logDet = True
profile = False
treinos = 10 #1200

# =================================================


### temporal pooler classes
#========================
class Segment:
    ''' In the temporal pooler, each segment can have its own (possibly large) list of
        potential synapses.
        To avoid excessive memory usage, we randomly add active synapses to each segment
        during learning.
    '''
    def __init__(self):
        # {column: permanence}
        self.potentialSynapses = {}

        ''' indicates whether the segment predicts feed-forward input on the next time step '''
        self.sequenceSegment = 0 


    def segmentActiveActiveState(self):
        ''' segmentActive(s, t, state): This routine returns true if the number of connected
            synapses on segment s that are active due to the active state at time t is greater
            than activationThreshold.
        '''
        actst = [col.id for col in self.potentialSynapses if self.potentialSynapses[col] > connectedPerm and col.overlap >= minOverlap]

        if logDet:
            print "segmentActiveActiveState():", actst

        qtact = len(actst)
        if qtact <= activationThreshold:
            qtact = 0
        return qtact
                  

#========================
class Cell:
    ''' Each cell is indexed using two numbers: a column index, c, and a cell index, i.
        Cells maintain a list of dendrite segments, where each segment contains a list of
        synapses plus a permanence value for each synapse. Changes to a cell's synapses
        are marked as temporary until the cell becomes active from feed-forward input.
        These temporary changes are maintained in segmentUpdateList.
        The pseudocode also uses a small state machine to keep track of the cell states at
        different time steps. We maintain three different states for each cell.
        The arrays activeState and predictiveState keep track of the active and predictive
        states of each cell at each time step. The array learnState determines which cell
        outputs are used during learning. When an input is unexpected, all the cells in a
        particular column become active in the same time step. Only one of these cells (the cell
        that best matches the input) has its learnState turned on. We only add synapses from
        cells that have learnState set to one (this avoids overrepresenting a fully active
        column in dendritic segments).

        Each distal dendrite segment also has an associated set of potential synapses.
        The set of potential synapses is a subset of all the cells in a region.
        As the segment learns, it increases or decreases the permanence value of all
        its potential synapses. Only those potential synapses that are above a threshold
        are valid.
    '''
    def __init__(self):
        # lista de Segment()
        self.dendriteSegments = []
        self.segmentUpdateList = []
        self.activeState = 0
        self.predictiveState = 0
        self.learnState = 0

    def getActiveSegmentActiveState(self):
        ''' for a given cell, return a segment index such that segmentActive(s,t, state) is true.
            If multiple segments are active, sequence segments are given preference.
            Otherwise, segments with most activity are given preference.
        '''
        acts = [(s.segmentActiveActiveState(), s.sequenceSegment)  for s in self.dendriteSegments]
        if logDet: print "getActiveSegmentActiveState():", acts
        if len(acts) > 0:
            return acts[0]
        else:
            return None
            
#========================
class TemporalPooler:
    def __init__(self):
        self.none = None


    def phase_1(self, actc):
        ''' The first phase calculates the active state for each cell.
            For each winning column we determine which cells should become active.
            If the bottom-up input was predicted by any cell (i.e. its predictiveState
            was 1 due to a sequence segment in the previous time step), then those cells
            become active. If the bottom-up input was unexpected (i.e. no cells
            had predictiveState output on), then each cell in the column becomes active.

            \param actc: list of active columns
        '''
        buPredicted = False
        lcChosen = False
        
        for column in actc:
            for cell in column.cells:
                if cell.predictiveState != 0:
                    s = cell.getActiveSegmentActiveState()
                    if s.sequenceSegment:
                        buPredicted = True
                        cell.activeState = 1
                    
            

###### spatial pooler classes
#========================
class History:
    def __init__(self):
        self.activeHistory = []  # for N past iterations, when column is active after inhibition
        self.overlapHistory = [] # for N past iterations, when overlap > minOverlap


#========================
class Column:    
    proxId = 0
    def __init__(self):
        self.id = Column.proxId
        Column.proxId += 1

        #information on proximal dendrite
        # {input: permanence}
        self.potentialSynapses = {}
        self.overlap = 0
        self.boost = 1.0
        self.activeHistory = []  # for N past iterations, when column is active after inhibition
        self.overlapHistory = [] # for N past iterations, when overlap > minOverlap
        self.activeDutyCycle = 1.0
        self.overlapDutyCycle = 1.0

        # information on cells and distal dendrites
        self.cells = [Cell() for i in range(cellsPerColumn)]

    def dump(self):
        print "potentialSynapses:"
        for ps in self.potentialSynapses:
            print "  ", ps, self.potentialSynapses[ps]
        print "overlap:", self.overlap
        print "boost:", self.boost
            

    def __repr__(self):
        ps = [(p, '%.2f'%self.potentialSynapses[p]) for p in self.potentialSynapses]
        return 'ovr:%d boost:%.1f syn:%s\nactiveHist:%s\noverlapHistory:%s'%(self.overlap, self.boost, str(ps), str(self.activeHistory), str(self.overlapHistory))

    def __str__(self_):
        self.__repr__()



# -------------------------------------------------
def limit(val, minv, maxv):
    return min(max(val,minv),maxv)
        
# =================================================
class SpatialPooler:
    def __init__(self):
        self.none = None
        
    # -------------------------------------------------
    def randomCenter(self):
        dx = 1.0-math.sin(random.random()*math.pi/2)        
        if random.random() < 0.5:
            dx = -dx

        return dx


    # -------------------------------------------------
    def initialization(self, dperm=0.1, npsyn=5, radius=5):
        ''' Prior to receiving any inputs, the region is initialized by computing a
        list of initial potential synapses for each column. This consists of a random set
        of inputs selected from the input space. Each input is represented by a synapse
        and assigned a random permanence value. The random permanence values are chosen
        with two criteria.
        * First, the values are chosen to be in a small range around
        connectedPerm (the minimum permanence value at which a synapse is considered "connected").
        This enables potential synapses to become connected (or disconnected) after a small number
        of training iterations.
        * Second, each column has a natural center over the input region, and
        the permanence values have a bias towards this center (they have higher values near the center).
        
        \param dperm: range around connectedPerm
        \param npsyn: number of potential synapses per column
        \param radius: size of region around natural center 
        '''
        sinpl = sinp[0]*sinp[1]
        scoll = scol[0]*scol[1]
        columns = []

        print "create %d columns with %d potentialSynapses for a space of %d inputs"%(scoll, npsyn, sinpl)
        
        for ic in range(scoll):
            col = Column()

            # center in input region
            cx = random.random()*sinp[0]
            cy = random.random()*sinp[1]

            if logDet: rep = 0
            while len(col.potentialSynapses) < npsyn:
                x = limit(int(cx + radius*self.randomCenter()), 0, sinp[0]-1)
                y = limit(int(cy + radius*self.randomCenter()), 0, sinp[1]-1)
                iinp = y*sinp[0]+x
                perm = connectedPerm - dperm + 2*dperm*random.random()
                if iinp in col.potentialSynapses:
                    if logDet: rep+=1
                col.potentialSynapses[iinp] = perm

            if logDet: print "col %d rep %d"%(ic,rep)
            columns.append(col)

        return columns
            

    # -------------------------------------------------
    def overlap(self, inputs, columns, minOvr = minOverlap):
        ''' Given an input vector, the first phase calculates the overlap of each column
            with that vector.
            The overlap for each column is simply the number of connected synapses with active inputs,
            multiplied by its boost. If this value is below minOverlap, we set the overlap score to
            zero.
        '''
        print "overlap %d inputs with %d columns"%(len(inputs),len(columns))
        for c in columns:
            ovr = 0
            for iinp in c.potentialSynapses:
                if c.potentialSynapses[iinp] > connectedPerm:
                    if inputs[iinp] > 0:
                        ovr += 1
            #if logDet: print "overlap column %d: %d"%(c.id,ovr)
            if ovr < minOvr:
                c.overlap = 0
                c.overlapHistory.append(0)
            else:
                c.overlap = ovr * c.boost
                c.overlapHistory.append(1)

            if len(c.overlapHistory) > historyLen:
                del c.overlapHistory[0]



    # -------------------------------------------------
    def inhibition(self, columns):
        ''' The second phase calculates which columns remain as winners after the inhibition step.
        desiredLocalActivity is a parameter that controls the number of columns that end up winning.
        For example, if desiredLocalActivity is 10, a column will be a winner if its overlap score is
        greater than the score of the 10'th highest column within its inhibition radius
        '''
        print "inhibition selects %d more active columns"%desiredLocalActivity
        colact = sorted([(c.overlap, c) for c in columns if c.overlap > 0], reverse=True)

        activeCols = colact[:desiredLocalActivity]

        if len(activeCols) > 0:
            if activeCols[-1][0] == 0:
                print "ERROR: activity 0 for %dth column"%desiredLocalActivity
        else:
            print "ERROR: 0 active columns!"

        actC = [c[1] for c in activeCols]
        for c in columns:
            if c in actC:
                c.activeHistory.append(1)
            else:
                c.activeHistory.append(0)
                
            if len(c.activeHistory) > historyLen:
                del c.activeHistory[0]
        

        return actC


    # -------------------------------------------------
    def learning1(self, actc, inputs):
        ''' This is part 1 of learning: updates of the permanence values of all synapses as
            necessary
            For winning columns, if a synapse is active, its permanence value is incremented, otherwise
            it is decremented. Permanence values are constrained to be between 0 and 1.
        '''

        print "learning1 %d active columns over %d inputs"%(len(actc),len(inputs))

        for column in actc:
            for iinp in column.potentialSynapses:
                if inputs[iinp] > 0:
                    column.potentialSynapses[iinp] = min(column.potentialSynapses[iinp]+permanenceInc, 1.0) 
                else:
                    column.potentialSynapses[iinp] = max(column.potentialSynapses[iinp]-permanenceDec, 0.0)


    # -------------------------------------------------
    def boostFunction(self, presentDC, minDC):
        ''' Returns the boost value of a column. The boost value is a scalar >= 1.
            If activeDutyCyle(c) is above minDutyCycle(c), the boost value is 1.
            The boost increases linearly once the column's activeDutyCyle starts
                 falling below its minDutyCycle.
        '''
        if presentDC >= minDC:
            return 1.0
        else:
            return 1.0+float(maxBoost-1.0)/minDC*(minDC-presentDC)
            

    # -------------------------------------------------
    def learning2(self, columns):
        ''' This is part 2 of learning: update of the boost and inhibition radius
            There are two separate boosting mechanisms in place to help a column learn connections.
            If a column does not win often enough (as measured by activeDutyCycle), its overall boost
            value is increased (line 30-32). Alternatively, if a column's connected synapses do not
            overlap well with any inputs often enough (as measured by overlapDutyCycle),
            its permanence values are boosted (line 34-36). Note: once learning is turned off,
            boost(c) is frozen.
            Finally, at the end of Phase 3 the inhibition radius is recomputed.
        '''
        lhis = len(columns[0].activeHistory)
        if lhis < historyLen:
            print "learning2: history = %d too short"%lhis
            # not enough history for boosting
            return

        # already has enough history for boosting

        # 1. update duty cycles
        maxDutyCycle = 0
        for c in columns:
            c.activeDutyCycle = sum(c.activeHistory)
            c.overlapDutyCycle = sum(c.overlapHistory)
            if maxDutyCycle < c.activeDutyCycle:
                maxDutyCycle = c.activeDutyCycle

        minDutyCycle = minDutyFactor * maxDutyCycle


        # 2. perform boosting
        boosted1 = 0
        boosted2 = 0

        dperm = 0.1*connectedPerm
        
        for c in columns:
            c.boost = self.boostFunction(c.activeDutyCycle, minDutyCycle)
            if c.boost > 1.001:
                boosted1 += 1

            if c.overlapDutyCycle < minDutyCycle:
    ##            if minDutyCycle == 0:
    ##                print "ERRO, c.overlapDutyCycle=",c.overlapDutyCycle
                boosted2 += 1
                for ps in c.potentialSynapses:
                    c.potentialSynapses[ps] = min(c.potentialSynapses[ps] + dperm, 1.0)
                    
            
        print "learning2 with maxDutyCycle = %d, minDutyCycle=%d, boosted=(%d,%d)"%(maxDutyCycle, int(minDutyCycle), boosted1, boosted2)
            

    # -------------------------------------------------
    def updateInputs(self, inputs, columns, learn=False):
        self.overlap(inputs, columns)
        actc = self.inhibition(columns)

        if learn:
            self.learning1(actc, inputs)
            self.learning2(columns)

        return actc


        
# -------------------------------------------------
def main():
    spatialPooler = None
    temporalPooler = None
    
    spatialPooler = SpatialPooler()
    temporalPooler = TemporalPooler()
    
    print datetime.datetime.now()
    if train:
        columns = spatialPooler.initialization(dperm, npsyn, radius)
        pickle.dump(columns, open('htm_cols.pick','wb'))
        
    else:
        print "reading column file"
        columns = pickle.load(open('htm_cols_t2.pick','rb'))
                    

    line=0
    line1 = [1 for i in range(sinp[0])]
    linetot = [0 for i in range(sinp[0]*sinp[1])]

    print datetime.datetime.now()
    if train:
        print "training"
        for i in range(treinos):
            inputs = linetot[0:line*sinp[0]]
            inputs.extend(line1)
            l2 = linetot[(line+1)*sinp[0]:]
            inputs.extend(l2)
            
            print "line=%d (%d)"%(line,len(inputs))
            line = (line+1)%sinp[1]

            actc = spatialPooler.updateInputs(inputs, columns, learn=True)

        print "dump columns"

##        self.activeHistory = []  # for N past iterations, when column is active after inhibition
##        self.overlapHistory = [] # for N past iterations, when overlap > minOverlap

        colHistory = []
        for c in columns:
            ch = History()
            ch.activeHistory = c.activeHistory
            ch.overlapHistory = c.overlapHistory
            colHistory.append(ch)

            c.activeHistory = []
            c.overlapHistory = []
        
        pickle.dump(columns, open('htm_cols_t2.pick','wb'))
        pickle.dump(colHistory, open('htm_cols_h2.pick','wb'))


    print datetime.datetime.now()
    if 1:
        print "testing"

        for i in range(20):
            inputs = linetot[0:line*sinp[0]]
            inputs.extend(line1)
            iinv = int(len(line1)*random.random())
            inputs[-iinv] = 0
            l2 = linetot[(line+1)*sinp[0]:]
            inputs.extend(l2)
            
            print "line=%d (%d)"%(line,len(inputs))
            line = (line+1)%sinp[1]
            

            for x in range(0,sinp[0]*sinp[1],sinp[0]):
                print inputs[x:x+sinp[0]]

            actc = spatialPooler.updateInputs(inputs, columns, learn=False)

            if temporalPooler != None:
                if logDet: print "calling tp phase 1"
                temporalPooler.phase_1(actc)

            cols = [(int(c.overlap), columns.index(c)) for c in actc]
            print sorted(cols, key=lambda x:x[1])
            print
            
        print


if __name__ == '__main__':
    print "INICIALIZANDO"
    print "train=",train
    print "logFile=",logFile
    print "logDet=",logDet
    print "profile=",profile
    print "treinos=",treinos

    a=raw_input("ENTER PARA CONTINUAR")

    bkpst = None

    if logFile:
        bkpst = sys.stdout
        sys.stdout = open("htm.log","at")


##    if profile:
##        profile.run('main()')
##    else:
##        main()

    main()
    
    if bkpst != None:
        sys.stdout.close()
        sys.stdout = bkpst
    
    a=raw_input("ENTER para encerrar")

