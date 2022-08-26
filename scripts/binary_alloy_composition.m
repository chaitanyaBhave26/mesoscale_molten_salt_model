clear;clc;

%%Ni-5Cr
m_cr = 0.05%20/100; % 5wt% Cr alloy

mol_ni = (1-m_cr)/58.693;
mol_cr = m_cr/51.996;

x_cr = mol_cr/(mol_ni+mol_cr);

M = x_cr*51.996 + (1-x_cr)*58.693; %g/mol, molar mass of approximated alloy
D = 8.7; %g/cm3

Va = (M/D)*(1e12/6.02214076e23) %um3/atom

%%Ni-20Cr

m_cr = 0.2%20/100; % 5wt% Cr alloy
m_c  = 0.05/100;

mol_ni = (1-m_cr)/58.693;
mol_cr = m_cr/51.996;

x_cr = mol_cr/(mol_ni+mol_cr);

M = x_cr*51.996 + (1-x_cr)*58.693; %g/mol, molar mass of approximated alloy
D = 8.57; %g/cm3

Va = (M/D)*(1e12/6.02214076e23) %um3/atom