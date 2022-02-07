# -*- coding: utf-8 -*-
"""
Streamlit app running locally to perform simulations with a RED model (co-flow,
NaCl).

Scientific literature on the model:
- Veerman et al. (2011), Chemical Engineering Journal
- Tedesco et al. (2015), Chemical Engineering Research and Design
- Simões et al. (2020), Desalination

D.Pintossi
2022-02-04
"""
import matplotlib.pyplot as plt
import streamlit as st
from model_mechanics.model_class import *

st.set_page_config(
    page_title='RED model (co-flow, NaCl)',
    layout='wide',
    initial_sidebar_state='expanded',
)


st.title('Reverse Electrodialysis process model')
st.markdown("*D.Pintossi, February '22*")

with st.sidebar:
    # Input residence time
    st.markdown('**Residence time**')
    single_residence_time = st.checkbox(
        'Equal residence time for seawater and river water',
        value=True,
        help='Check the box to have the same residence time for SW and RW',
    )
    if single_residence_time:
        tSW = st.slider(
            label='Residence time [s]:',
            min_value=10,
            max_value=100,
            step=2,
            help='Slide to set the residence time',
            key='single-residence-time',
        )
        tRW = tSW
    else:
        tSW = st.slider(
            label='Residence time for seawater [s]:',
            min_value=10,
            value=22.0,
            max_value=100,
            step=2,
            help='Slide to set the residence time',
            key='single-residence-time',
        )
        tRW = st.slider(
            label='Residence time for river water [s]:',
            min_value=10,
            value=22.0,
            max_value=100,
            step=2,
            help='Slide to set the residence time',
            key='single-residence-time',
        )

    # Geometry
    st.markdown('___')
    st.markdown('**Geometry**')
    L = st.slider(
        label='Length of the flow path [cm]:',
        min_value=6.5,
        value=22.0,
        step=0.5,
        max_value=44.0,
        help='Insert the length in centimeters',
        key='length',
    ) / 100
    W = st.slider(
        label='Width of the flow path [cm]:',
        min_value=6.5,
        value=22.0,
        step=0.5,
        max_value=44.0,
        help='Insert the width in centimeters',
        key='width',
    ) / 100
    d = st.slider(
        label='Thickness of the flow channel [μm]:',
        min_value=15,
        value=100,
        step=1,
        max_value=1000,
        help='Insert the channel thickness in micrometers',
        key='thickness',
    ) / 1000000

    # Constants
    Rblank = 0
    obstr = 1

    # Membranes
    st.markdown('___')
    st.markdown('**Membranes**')
    membranes = st.radio(
        label='Choose the membranes',
        options=['Ideal CEM/AEM', 'Fujifilm CEM/AEM type 10'],
        key='membrane-choice'
    )

    # Waters
    st.markdown('___')
    st.markdown('**Waters (BCs)**')
    cRW0 = st.number_input(
        label='Inlet NaCl concentration for seawater [mM]:',
        min_value=1.0,
        value=17.1,
        step=0.1,
        max_value=50.0,
        key='c0_sw',
    )
    cSW0 = st.number_input(
        label='Inlet NaCl concentration for river water [-]:',
        min_value=350,
        value=512,
        step=2,
        max_value=750,
        key='c0_rw',
    )

    # External load voltage
    st.markdown('___')
    st.markdown('**External load voltage**')
    Uload = st.slider(
        label='External load voltage per cell pair [mV]:',
        min_value=0,
        max_value=160,
        step=5,
        help='Slide to set the external load voltage per cell pair',
        key='Uload',
    ) / 1000

    # Number of intervals and iterations
    st.markdown('___')
    st.markdown('**Discretization parameters (BCs)**')
    intervals = st.slider(
        label='Number of intervals along each axis [-]:',
        min_value=20,
        value=100,
        step=20,
        max_value=1000,
        key='intervals',
    )

RED_model = CoFlow(
    tRW, tSW,
    W, L, d,
    membranes,
    Uload,
    cRW0, cSW0,
    intervals,
)

RED_model.solve_system()

st.markdown('___')
columns = st.columns(4)
with columns[0]:
    st.metric(
        label='Power',
        value=f'{RED_model.power:4.3f} W',
    )
with columns[1]:
    st.metric(
        label='Power density',
        value=f'{RED_model.pd_avg:4.3f} W/m2',
    )
with columns[2]:
    st.metric(
        label='Energy efficiency',
        value=f'{RED_model.eta_energy:4.1f} %',
    )
with columns[3]:
    st.metric(
        label='Thermodynamic efficiency',
        value=f'{RED_model.eta_thermo:4.1f} %',
    )
st.markdown('___')

# Plot concentrations
conc_fig, ax = plt.subplots(figsize=(6, 3))
ax.plot(RED_model.x, RED_model.cSW, label='[NaCl] in SW', color='#f63366')
ax.plot(RED_model.x, RED_model.cRW, label='[NaCl] in RW', color='#262730')
ax.set_xlabel('L [m]', fontweight='bold')
ax.set_ylabel('C [mM]', fontweight='bold')
ax.legend()

# Plot emf
emf_fig, ax2 = plt.subplots(figsize=(6, 3))
ax2.plot(RED_model.x, RED_model.E, label='emf', color='#f63366')
ax2.set_xlabel('L [m]', fontweight='bold')
ax2.set_ylabel('emf [V]', fontweight='bold')
ax2.set_ylim([0, 0.160])
ax2.legend()

columns2 = st.columns(3)
with columns2[1]:
    st.pyplot(fig=conc_fig)
    st.write('')
    st.pyplot(fig=emf_fig)
