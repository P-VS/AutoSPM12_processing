function f= acid_hysco_fermi(x,x0,kt)
% This is a fermi function that goes from 0 to 2 at the point x0. kt
% provides the steepness with which the transition is done.
% S.Mohammadi 2.10.2019

f= 2./(exp((abs(x)-x0).*kt)+1);