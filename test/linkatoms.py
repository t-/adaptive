from gmx_pdb import gmx_pdb
from gmx_top import gmx_topology
from math import *
from copy import *
import os

class read_pdb ():
    def __init__(self,pdb='mol.pdb',top='topol.top'):
        os.system('mkdir res') #qm jobs will go in here ;-)
        #--lamda stuff
        self.R1=1.0
        self.R2=6.0
        #--end
        self.pdb        = gmx_pdb(pdb)
        self.top        = gmx_topology()
        self.top.read_top(top)
        self.assign_bonds() #each atom gets a list of other atoms it is bound to
        #print self.pdb.atoms[1002].bonds
        self.qm_connect = self.assign_qmgroups()  # list of [ [groupmembers],[partners] ] with [partners]=[[groupmember],[outsideatom]] so [ [1,2,3], [[2,5],[3,7]]
        #print self.pdb.atoms[1002].qm_group
        #print self.qm_connect[self.pdb.atoms[1002].qm_group]
        self.center          = [self.pdb.atoms[1008-1].x,self.pdb.atoms[1008-1].y,self.pdb.atoms[1008-1].z]
        self.azone,self.zone_id_list,self.azone_center_id = self.find_azone(self.center,radius=self.R2,center_atom_number=1008)
        self.zones = self.create_adaptive_zones()
        self.mm=[]
        self.mm_ids=[]
        self.qm=[]
        self.qm_ids=[]
        self.input_overhead=[]
        self.read_gin()
        self.create_zones_top()
        self.fort7=self.read_fort7()
        
    def assign_bonds(self):
        for at in self.top.bonds_list:
            try:
                i =(at.ai)
                ip=(at.ai-1)
                j =(at.aj)
                self.pdb.atoms[ip].bonds[i].append(j)
            except Exception:
                self.pdb.atoms[ip].bonds[i] = [j]
            try:
                j =(at.aj)
                jp=(at.aj-1)
                i =(at.ai)
                self.pdb.atoms[jp].bonds[j].append(i)
            except Exception:
                self.pdb.atoms[jp].bonds[j] = [i]
        #for i in self.pdb.residues:
        #    residue_types[i[0].res]=''
        #print residue_types.keys()
        
    def assign_qmgroups(self):
        resgroup=[]
        for res in xrange(len(self.pdb.residues)):
            tmpgroup=[]
            for atm in self.pdb.residues[res]:
                if atm.name not in ['N','H','C','O','CA','HA','LA'] and atm.qm_group==-1:#'OW','HW1','HW2'
                    tmpgroup.append(atm.nr)
                #    if atm.name not in ['OW','HW1','HW2']:print atm.nr,atm.name #<-- here, everything still works ;-)
            for i in tmpgroup:
                self.pdb.atoms[i-1].qm_group=len(resgroup)
            if tmpgroup:
                resgroup.append(tmpgroup)
            
        qmconnect=[]
        for group in resgroup:
            partners = []
            partners_bonds =[]
            for atm in group:
                try:
                    for partner in self.pdb.atoms[atm-1].bonds[atm]:
                        if partner not in group:
                            partners.append(partner) #####################
                except Exception:a=0;
            qmconnect.append([group,partners])
        return qmconnect
    
    def find_azone(self,center,radius,center_atom_number):
        #creates and completes a zone of radius r around a center atom
        zone_id_list={}
        azone_center_id=-1
        c=center
        azone={}
        for i in self.pdb.atoms:    #cuts QM_zone
            if (i.x-c[0])**2 + (i.y-c[1])**2 + (i.z-c[2])**2 <= radius**2:
                if i.qm_group!=-1:
                    azone[i.nr] = ''
        #-----writes_pdb_for_illustrational_purposes-----#
        f = open('azone_pre.pdb','w')
        tmp=[]
        for i in self.pdb.atoms:
            if i.nr in azone.keys():
                tmp.append(i)
        for i in tmp:
            print >>f, i.cprint()
        f.close()
        #-----END-----#

        #here, the QM_zone is completed (or T-Zone in the adaptive case)
        for key in azone.keys():        #loops over QM_zone atom numbers
            atm = self.pdb.atoms[key-1] #call pdb atom with matching atom number
            if atm.qm_group!= -1 :      #exclude all backbone atoms (also all PRO & GLY residues)
                zone_id_list[atm.qm_group]='' #creates a list of qm_group ids which are in the A&T zone
                k = self.qm_connect[atm.qm_group][0]
                for atom in k:
                        if atom==center_atom_number:
                            azone_center_id=atm.qm_group    #gets the qm_group id of the A-zone qm_group
                        azone[atom]=''
        
        #-----writes_pdb_for_illustrational_purposes-----#
        f = open('azone_post.pdb','w')
        tmp=[]
        
        for atm in azone.keys():
            tmp.append(self.pdb.atoms[atm-1])
        for i in tmp:
            print >>f, i.cprint()
        f.close()
        #-----END-----#
        return azone,zone_id_list,azone_center_id
    
    def calc_lamda(self,r):
        if r<self.R1:
            return 0.0
        elif r<=self.R2 and r>=self.R1:
            return (r-self.R1)**2.0 * (3*self.R2 - self.R1 - 2*r) / (self.R2 - self.R1)**3.0
        elif r>self.R2:
            return 1.0

    #old
    def create_adaptive_zones(self):
        lamda_list = {}
        for zone in self.zone_id_list:
            com=[0.0,0.0,0.0]
            count = 0.0 
            for atomnr in self.qm_connect[zone][0]:
                if self.pdb.atoms[atomnr-1].name.find('H')<=0:
                    com[0]=com[0]+self.pdb.atoms[atomnr-1].x
                    com[1]=com[1]+self.pdb.atoms[atomnr-1].y
                    com[2]=com[2]+self.pdb.atoms[atomnr-1].z
                    count=count+1.0
            if count>0.0:
                com[0]=com[0]/count
                com[1]=com[1]/count
                com[2]=com[2]/count
            if zone!=self.azone_center_id:
                lamda = self.calc_lamda( sqrt( (com[0]-self.center[0])**2 + (com[1]-self.center[1])**2 + (com[2]-self.center[2])**2  ) )
                #print '%.3f' %lamda, '%.2f' %sqrt( (com[0]-self.center[0])**2 + (com[1]-self.center[1])**2 + (com[2]-self.center[2])**2  )
            elif zone==self.azone_center_id:
                lamda = 0.0
            self.qm_connect[zone].append(com)
            self.qm_connect[zone].append(lamda)
            lamda_list[lamda]=zone
        lamdas= lamda_list.keys()
        lamdas.sort()
        #print lamdas
        zones=[]
        for i in xrange(len(lamdas)):
            tmp=[]
            if lamdas[i]<1:
                sigma=lamdas[i+1]-lamdas[i]
                for part in lamdas[:i+1]:
                    tmp.append(lamda_list[part])
                zones.append([tmp,sigma])
        #print zones
        return zones
    
    #old
    def create_zones_top_working(self):
        for zone in self.zones: #get all atoms in all zones
            top = deepcopy(self.top)                        ######## THIS IS EXTREMELY SLOW... here one should us an 'undo' function 
            pdb = deepcopy(self.pdb)                        ######## THIS IS EXTREMELY SLOW... here one should us an 'undo' function
            print 'creating zone',zone[0],'of sigma=',zone[1]
            qm_zone_atoms = []
            for zone_number in zone[0]:
                for zone_atom in self.qm_connect[zone_number][0]:
                    qm_zone_atoms.append(zone_atom)

            group_neighbours={} #get neighbours of all zones (e.i. atoms bound to the zone but not member of the zone)
#            for zone in self.zones:
            for zone_number in zone[0]:
                for zone_neighbour in self.qm_connect[zone_number][1]:  
                    group_neighbours[zone_neighbour] = ''

            links=[]
            for i in group_neighbours.keys(): # reduce to only neighbours of zones that are _not_ part of another included zone!
                for j in pdb.atoms[i-1].bonds[i]:
                    if j in qm_zone_atoms:
                        print 'placing link atoms on bonds:', j,'-',i
                        links.append([j,i]) #link[qm_atom,mm_atom]
                        
            #==TOPology==# bonds section fix
            for bond in top.bonds_list: #set all the bonds in the qm zone to function type 5
                if bond.ai in qm_zone_atoms: 
                    bond.funct = 5
                if bond.aj in qm_zone_atoms:
                    bond.funct = 5
                for link in links: # remove explicit bonds between qm and mm zone
                    if (bond.ai==link[0] and bond.aj==link[1]) or (bond.ai==link[1] and bond.aj==link[0]):
                        bond.funct = -1
            #==TOPology==# atoms section fix & vsites section fix & constraints section fix
            top.virtual_sites2_meta.append('[ virtualsites2 ]')             
            top.constraints_meta.append('[ constraints ]')                  
            for i in links:
            #pdb_file
                tmp4 = pdb.pdb_atom()
                tmp4.nr=str(len(top.atoms_list)+1)
                tmp4.name='LA'
                tmp4.res='XXX'
                tmp4.seq=str(top.atoms_list[-1].resnr+1)
                tmp4.x='0.0'
                tmp4.y='0.0'
                tmp4.z='0.0'
                tmp4.occup='1.0'
                tmp4.tempfac='0.0'
                pdb.atoms.insert(len(top.atoms_list),tmp4)
            #linkatoms
                tmp = top.atoms(\
                               nr       = str(len(top.atoms_list)+1) ,\
                               type     = 'LA',\
                               resnr    = str(top.atoms_list[-1].resnr+1),\
                               residue  = 'XXX',\
                               atom     = 'LA',\
                               cgnr     = str(top.atoms_list[-1].cgnr+1),\
                               charge   = '0.0',\
                               mass     = '0.0' )
                top.atoms_list.append(tmp)
            #vsites
                tmp2 = top.virtual_sites2(\
                                site    = str(len(top.atoms_list)+1),\
                                qmatom  = i[0],\
                                mmatom  = i[1],\
                                type    = 1,\
                                a       = 0.65)
                top.virtual_sites2_list.append(tmp2)
            #constraints
                tmp3 =top.constraints(\
                                qmatom  = i[0],\
                                mmatom  = i[1],\
                                type    = 2,\
                                a       = 0.153)
                top.constraints_list.append(tmp3)
            print 'writing:','link_'+str(zone[1])+'.top/.pdb'
            top.write_top('link_'+str(zone[1])+'.top')
            pdb.write_pdb('link_'+str(zone[1])+'.pdb')
        
    def create_zones_top(self):
        for zone in self.zones: #get all atoms in all zones
            tmp_zone_id_list={}
            print 'creating zone',zone[0],'of sigma=',zone[1]
            qm_zone_atoms = []
            for zone_number in zone[0]:
                for zone_atom in self.qm_connect[zone_number][0]:
                    qm_zone_atoms.append(zone_atom)

            group_neighbours={} #get neighbours of all zones (e.i. atoms bound to the zone but not member of the zone)
#            for zone in self.zones:
            for zone_number in zone[0]:
                for zone_neighbour in self.qm_connect[zone_number][1]:  
                    group_neighbours[zone_neighbour] = ''
            links=[]
            for i in group_neighbours.keys(): # reduce to only neighbours of zones that are _not_ part of another included zone!
                atm = self.pdb.atoms[i-1]
                for j in atm.bonds[i]:
                    #if 1000<j<1010:print i,j
                    if j in qm_zone_atoms:
                        #if 699<j<1010:print i,j
                        print 'placing link atoms on bonds:', j,'-',i
                        links.append([j,i,self.dummy(j,i)]) #link[qm_atom,mm_atom]
            g_qm_new=[]
            g_mm_new=[]
            qm_residues={}            
            qm_charge=0.0            
            for i in qm_zone_atoms: # create the lines for the QM region in the gaussian input file
                qm_residues[self.pdb.atoms[i-1].seq]=''
                if self.pdb.atoms[i-1].name.strip('123456789')[0] =='C':
                    protons=6
                if self.pdb.atoms[i-1].name.strip('123456789')[0] =='N':
                    protons=7
                if self.pdb.atoms[i-1].name.strip('123456789')[0] =='O':
                    protons=8
                if self.pdb.atoms[i-1].name.strip('123456789')[0] =='H':
                    protons=1
                if self.pdb.atoms[i-1].name.strip('123456789')[0] =='P':
                    protons=15
                if self.pdb.atoms[i-1].name.strip('123456789')[0] =='S':
                    protons=16
                try:
                    atm=self.pdb.atoms[i-1]
                    line='%3d  %3.7f  %4.7f  %4.7f'%(protons,atm.x,atm.y,atm.z)
                    g_qm_new.append(line)
                except Exception:
                    a=0
                tmp_zone_id_list[i]=''
            dum_list=[]        
            for i in links:
                line='  1  %3.7f  %4.7f  %4.7f'%(i[2][0],i[2][1],i[2][2])
                dum_list.append(len(g_qm_new))
                g_qm_new.append(line)
            zone.append(dum_list)
            for i in qm_residues.keys():    #calculate the charge of the new adaptive region for gaussian input file
                q=0.0
                for atom in self.pdb.residues[i-1]:
                    nr=atom.nr-1
                    atm=self.pdb.atoms[nr]
                    if atm.name not in ['OW','HW1','HW2','LA','NA','CL']:
                        q=q+self.top.atoms_list[nr].charge
                    if atm.name in ['OW']:
                        q=q+(-0.834)
                    if atm.name in ['HW1','HW2']:
                        q=q+(0.417)
                    if atm.name in ['NA']:
                        q=q+(1.0)
                    if atm.name in ['CL']:
                        q=q+(-1.0)
                qm_charge=qm_charge+q
            print 'charge of adaptive QM region:%1.0f'%qm_charge

            for i in self.mm_ids:
                if i not in tmp_zone_id_list:
                    atm=self.pdb.atoms[i-1]
                    nr=atm.nr-1
                    if atm.name not in ['OW','HW1','HW2','LA','NA','CL']:
                        charge=self.top.atoms_list[nr].charge
                    if atm.name in ['OW']:
                        charge=-0.834
                    if atm.name in ['HW1','HW2']:
                        charge=0.417
                    if atm.name in ['NA']:
                        charge=1.0
                    if atm.name in ['CL']:
                        charge=-1.0
                    line='%.7f  %4.7f  %4.7f  %4.4f'%(atm.x,atm.y,atm.z,charge)
                    g_mm_new.append(line)

            #print gaussian input file!
            os.system('mkdir res/'+str(zone[1]))
            os.system('cp run_g03.sh res/'+str(zone[1])+'/')
            f=open('res/'+str(zone[1])+'/mol.com','w')
            for i in self.input_overhead:
                print >>f, i
            print >>f, ''
            print >>f, '%1.0f 1'%qm_charge
            for i in g_qm_new:
                print >>f, i
            print >>f,''
            for i in g_mm_new:
                print >>f,i
            print >>f, ''
            f.close()
            os.system('cd res/'+str(zone[1])+'/; sh run_g03.sh')
            
    def dummy(self,qm,mm):
        qatm=self.pdb.atoms[qm-1]
        r_qm=[qatm.x,qatm.y,qatm.z]
        matm=self.pdb.atoms[mm-1]
        r_mm=[matm.x,matm.y,matm.z]
        r_c=[0.0,0.0,0.0]
        #rmm-rqm
        r_c=[r_mm[0]-r_qm[0],r_mm[1]-r_qm[1],r_mm[2]-r_qm[2]]
        r_c_abs=sqrt( r_c[0]**2 + r_c[1]**2 + r_c[2]**2)
        r_c=[r_c[0]*1.1/r_c_abs,r_c[1]*1.1/r_c_abs,r_c[2]*1.1/r_c_abs]
        dum = [ r_qm[0]+r_c[0],r_qm[1]+r_c[1], r_qm[2]+r_c[2] ]
        return dum
                                    
    def read_fort7(self):
            fort=[]
            for zone in self.zones: #get all atoms in all zones
                tmp=[]
                new=0
                energy=0
                for line in open('res/'+str(zone[1])+'/fort.7','r'):
                    line=line.strip('\n').split(' ')
                    for i in range(line.count('')):line.remove('')
                    if len(line)==1 and new==0:
                        new=1
                        energy = float(line[0])
                    if len(line)==3:
                        tmp.append([float(line[0])*zone[1],float(line[1])*zone[1],float(line[2])*zone[1]]) #scale all the forces by sigma!! 
                fort.append([tmp,energy,zone[2]])
                
            new=[]
            for anr in xrange(len(fort[0][0])): #loop over atoms 
                for znr in xrange(len(fort)):#loop over force on this atom from different zones
                    nr = anr
                    #fort[zone][linkatompositions][firstlinkatom]
                    try:
                        if anr>=fort[znr][2][1]:        #there is one link atom that has been present on the PYP chromophore from the beginning fort[znr][2][0]. This, we keep! all the other ones (fort[znr][2][1,2,..,etc.]) get dropped
                            nr=anr+len(fort[znr][2])-1  #len(fort[znr][2])-1 = number of linkatoms excluding the initial one. shift the index numbers away from the link atom, however include the one linkatom from the A-region as also present in gromacs!
                    except Exception: 0; 
                    cur=fort[0][0][anr]
                    cur=[cur[0]+fort[znr][0][anr][0],cur[1]+fort[znr][0][anr][1],cur[2]+fort[znr][0][anr][2]]
                new.append(cur)
            f=open('fort.8','w')
            print >>f,'   -0000.000000000000   '
            for i in new:
                print >>f,'  %5.10e  %5.10e  %5.10e'%(i[0],i[1],i[2])
            print >>f,'   -1   '
            f.close()
                
                    
    def read_gin(self):
        blanks=0
        input_overhead=[]
        qm=[]
        mm=[]
        qm_ids=[]
        mm_ids={}
        qm_charge=0.0
        for line in open('input.com','r'):
            if line=='\n':
                blanks=blanks+1
            if blanks<2:
                line=line.replace('units=bohr','').strip('\n')
                input_overhead.append(line)
            if blanks==2:
                line=line.strip('\n')
                line=line.split(' ')
                for i in range(line.count('')): line.remove(''); #remove blanks from array
                if len(line)>2:
                    qm.append(line)
            if blanks==3:
                line=line.strip('\n')
                line=line.split(' ')
                for i in range(line.count('')): line.remove(''); #remove blanks from array
                if len(line)>2:
                        mm.append(line)
        for line in open('mm_ids','r'):
            mm_ids[int(line)]=''
        for line in open('qm_ids','r'):
            qm_ids.append(int(line))
        self.mm=mm
        self.mm_ids=mm_ids
        self.qm=qm
        self.qm_ids=qm_ids
        self.input_overhead=input_overhead
p =read_pdb()