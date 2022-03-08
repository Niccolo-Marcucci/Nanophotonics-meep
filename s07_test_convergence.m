% clc
clear all
close all

addpath('Near-to-Far-Field-Transformation')
t = ['10'; '20'; '30'; '40'];
name = 'polSplitter_6ec58eed79_design_TM_gd3_buriedDBR_onSiO2_positive_N9_Dphi60_sigma-1_charge0-res133_220114-192330_nearfield_t';

for i = 1:4

f1 = [name, (t(i,:)),'.0.mat'];

fields = load(['data/',f1]);

volume = [fields.Lx, fields.Ly, 0]*1e-6;

Ex1 = fields.Ex;
Ey1 = fields.Ey;
Ex1 = Ex1(2:end,:);
Ey1 = Ey1(:,2:end);

N = size(Ex1);
x = linspace(-1,1,N(1))*volume(1);
y = linspace(-1,1,N(2))*volume(2);
[X,Y] = meshgrid(y,x);
c0 = physconst('lightspeed');

wavelength = 0.57*1e-6;

freq = c0/wavelength;
obj_near = Field2D(X, Y, Ex1, Ey1, freq, 'm');
figure
obj_near.plotNearField()


f2 = [name, (t(i+1,:)),'.0.mat'];
fields = load(['data/',f2]);

Ex2 = fields.Ex;
Ey2 = fields.Ey;
Ex2 = Ex2(2:end,:);
Ey2 = Ey2(:,2:end);

E1 = sum(sum(abs(Ex1))) + sum(sum(abs(Ey1)));
E2 = sum(sum(abs(Ex2))) + sum(sum(abs(Ey2)));

DE = log(abs( E1 - E2 ))
relDE1 = log(abs( E1 - E2 ) / E1)
relDE2 = log(abs( E1 - E2 ) / E2)
end
