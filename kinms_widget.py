#!/usr/bin/env python3
# coding: utf-8
import numpy as np
from kinms import KinMS
# import matplotlib.pyplot as plt
# from gastimator import corner_plot
from kinms_fitter.sb_profs import sb_profs
from kinms_fitter.velocity_profs import velocity_profs
import matplotlib.pyplot as plt
from kinms.utils.KinMS_figures import KinMS_plotter
import streamlit as st
from astropy.io import fits

class model_maker:
    def __init__(self,velprof,sbprof):
        self.pa = 270
        #self.vsys
        self.inc = 60
        self.totflux = 1
        self.veldisp = 8
        self.vel_profile = velprof
        self.sb_profile = sbprof
        self.beam = 3.
        self.cellsize = self.beam/3.
        
        self.nx = 128
        self.ny = 128
        self.nv = 100
        self.dv = 10
        self.nSamps=5e5
        self.bunit = 'Jy/beam'
        self.sbRad = np.arange(1e-5,self.nx*self.cellsize*2,self.beam/6.)
        
    def make(self,fileName='',plot=True):
        pa=self.pa
        #xc=param[1]
        #yc=param[2]
        #vsys=self.vsys
        inc=self.inc
        totflux=self.totflux
        veldisp=self.veldisp
        #phasecen=[xc,yc]
        self.bmaj,self.bmin,self.bpa = self.beam, self.beam, 0
        self.n_sbvars = np.sum([i.freeparams for i in self.sb_profile])
        self.n_velvars = np.sum([i.freeparams for i in self.vel_profile])
        
        vrad=velocity_profs.eval(self.vel_profile,self.sbRad,np.concatenate([i.guess for i in self.vel_profile]).ravel(),inc=inc)
        sbprof=sb_profs.eval(self.sb_profile,self.sbRad,np.concatenate([i.guess for i in self.sb_profile]).ravel())
        # myargs={'phaseCent': phasecen}
        #st.write([self.nx*self.cellsize,self.ny*self.cellsize,self.nv*self.dv,self.cellsize,self.dv])
        
        cube=KinMS(self.nx*self.cellsize,self.ny*self.cellsize,self.nv*self.dv,self.cellsize,self.dv,\
                 [self.bmaj,self.bmin,self.bpa],fixSeed=True,nSamps=self.nSamps).model_cube(inc,sbProf=sbprof,sbRad=self.sbRad,velRad=self.sbRad,velProf=vrad,gasSigma=veldisp,\
                 intFlux=totflux,posAng=pa,fileName=fileName,bunit=self.bunit)
              
        #st.text([self.nx*self.cellsize,self.ny*self.cellsize,self.nv*self.dv,self.cellsize,self.dv,bem[self.bmaj,self.bmin,self.bpa]])        
        
        plotter=KinMS_plotter(cube,self.nx*self.cellsize,self.ny*self.cellsize,self.nv*self.dv,self.cellsize,self.dv,[self.bmaj,self.bmin,self.bpa], posang=pa)
        plotter.makeplots()

        return plotter 
        

st.set_page_config(
     page_title="KinMS Kinematics visualiser",
     layout="centered",
     initial_sidebar_state="expanded",
     # menu_items={
     #     'Get Help': 'https://www.extremelycoolapp.com/help',
     #     'Report a bug': "https://www.extremelycoolapp.com/bug",
     #     'About': "# This is a header. This is an *extremely* cool app!"
     # }
 )


# st.markdown(
# """
# <style>
# [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
# width: 30rem;
# }
# [data-testid="stSidebar"][aria-expanded="false"] > div:first-child {
# width: 30rem;
# }
# </style>
# """,
# unsafe_allow_html=True
# )        
             
#st.image("https://raw.githubusercontent.com/TimothyADavis/KinMSpy/master/kinms/docs/Logo.png",width=200) 
st.sidebar.title('KinMS kinematics visualiser')

#with st.sidebar.expander('Telescope Setup', expanded=False):
tab1, tab2 = st.sidebar.tabs(["Object Setup","Telescope Setup"])
with tab2:    
    nxarc=st.number_input('Cube size (arcsec)',value=64.,step=1.,key='nx',max_value=120.,min_value=10.)  
    #nyarc=st.sidebar.number_input('Y size (arcsec)',value=64.,step=1.,key='ny')
    nchans=st.number_input('Number of channels',value=100,step=1,key='nchans',max_value=200,min_value=5)  
    dv=st.number_input('Channel width (km/s)',value=10,step=5,key='dv',max_value=120,min_value=1)  
    beam=st.slider('Beamsize (arcsec)', min_value=0.3, max_value=10., value=3.)

with tab1:
    pa=st.slider('PA (deg)', min_value=0, max_value=360, value=270)
    inc=st.slider('Inclination (deg)', min_value=0, max_value=90, value=45)

    distance=st.number_input('Distance (Mpc)',value=16.5,step=0.5,key='dist')   
    cellsize=beam/3.




col1,col2 = st.columns(2)
with col1:
    my_slot1 = st.empty()
    with st.expander('Edit SB profile', expanded=False):

        method_list = np.array(['expdisk','gaussian','cutoff'])
        selected=st.multiselect('Model Components', method_list,default='expdisk')
        my_slot1.subheader('SB Model: '+' + '.join(selected))
        my_sb_model=[]
        default_fixed=None
        for model in selected:
            st.write(model+': Required Parameters')
            if model=='expdisk':
                if 'gaussian' in selected:
                    needpars=[0,1]
                    bestguess=[0.5,10.]
                    bestmin=[0.0,0.0]
                    bestmax=[1.0,nxarc]
                else:
                    needpars=[1.]
                    bestguess=[10.]
                    bestmin=[0.0]
                    bestmax=[nxarc]
            if model=='gaussian':
                if 'expdisk' in selected:
                    needpars=[0.,1.,2.]
                    bestguess=[0.5,5.,1.]
                    bestmin=[0.0,0.01,0.01]
                    bestmax=[1.0,nxarc*1.0,nxarc*0.5]
                else:
                    needpars=[1.,2.]
                    bestguess=[5.,1.]
                    bestmin=[0.0,0.01]
                    bestmax=[nxarc,nxarc*0.5]
            if model=='cutoff':
                needpars=[0.,1.]
                bestguess=[0.0,5.0]
                bestmin=[0.0,0.01]
                bestmax=[0.0,nxarc*0.5]
                default_fixed=[True,False]
            
            thismodel=getattr(sb_profs, model)(needpars,needpars,needpars)
            labels=thismodel.labels
            units=thismodel.units

            required_vals=[]
            min_vals=[]
            max_vals=[]
            fixed=[]
            for index in range(0,len(labels)):
                st.caption(labels[index]+' ('+units[index]+')')
                required_vals.append(st.number_input('Guess',key=labels[index],step=0.1,value=bestguess[index]))




            #text=st.sidebar.subheader('Optional Parameters')
            my_sb_model.append(getattr(sb_profs, model)(required_vals,min_vals,max_vals,fixed=fixed))

with col2:
    my_slot2 = st.empty()
    with st.expander('Edit velocity profile', expanded=False):
        method_list = np.array(['arctan','bulge_disc','keplarian','nfw','sersic'])
        selected2=st.multiselect('Model Components', method_list,default='arctan')
        my_slot2.subheader('Vcirc Model: '+' + '.join(selected2))
        my_vel_model=[]
        for model in selected2:
            if model == 'arctan':
                bestguess=[200.,1.]
                bestmin=[0.0,0.0,0.0,0.0]
                bestmax=[0.0,0.0,0.0,0.0]
            if model == 'bulge_disc':
                bestguess=[10.5,10.,5.,0.3]
                bestmin=[0.0,0.0,0.0,0.0]
                bestmax=[0.0,0.0,0.0,0.0]
            if model == 'keplarian':
                bestguess=[8.]
                bestmin=[0.0,0.0,0.0,0.0]
                bestmax=[0.0,0.0,0.0,0.0]
            if model == 'nfw':
                bestguess=[12.,100.]
                bestmin=[0.0,0.0,0.0,0.0]
                bestmax=[0.0,0.0,0.0,0.0]
            if model == 'sersic':
                bestguess=[10.5,10.,2.]
                bestmin=[0.0,0.0,0.0,0.0]
                bestmax=[0.0,0.0,0.0,0.0]


            st.write(model+': Required Parameters')
            try:
                thismodel=getattr(velocity_profs, model)(np.zeros(5),np.zeros(5),np.zeros(5))
            except:
                thismodel=getattr(velocity_profs, model)(distance,np.zeros(5),np.zeros(5),np.zeros(5))
            labels=thismodel.labels
            units=thismodel.units

            required_vals=[]
            min_vals=[]
            max_vals=[]
            for index in range(0,len(labels)):
                st.caption(labels[index]+' ('+units[index]+')')
                required_vals.append(st.number_input('Guess',key=labels[index],step=0.1,value=bestguess[index]))

            try:
                my_vel_model.append(getattr(velocity_profs, model)(required_vals,bestmin,bestmax))
            except:
                #breakpoint()
                my_vel_model.append(getattr(velocity_profs, model)(distance,guesses=required_vals,minimums=bestmin,maximums=bestmax))


plot_slot = st.empty()

if len(my_vel_model)>0 and len(my_sb_model)>0:

    model=model_maker(my_vel_model,my_sb_model)
    model.pa = pa
    model.inc = inc
    model.beam = beam
    model.cellsize = cellsize
    model.dv = dv
    model.nv= nchans
    model.nx= np.floor(nxarc/cellsize)
    model.ny= model.nx

    with st.spinner('Working...'):
        figmade=model.make()
    
        #buf = BytesIO()
        #figmade.figure.savefig(buf, format="png",bbox_inches='tight')
    
    
        #st.write(figmade)
        plot_slot.pyplot(figmade.figure)
        
        