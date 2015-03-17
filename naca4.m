% C***********************************************************************
% C    Module:  naca.f
% C
% C    Copyright (C) 2000 Mark Drela
% C
% C    This program is free software; you can redistribute it and/or modify
% C    it under the terms of the GNU General Public License as published by
% C    the Free Software Foundation; either version 2 of the License, or
% C    (at your option) any later version.
% C
% C    This program is distributed in the hope that it will be useful,
% C    but WITHOUT ANY WARRANTY; without even the implied warranty of
% C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% C    GNU General Public License for more details.
% C
% C    You should have received a copy of the GNU General Public License
% C    along with this program; if not, write to the Free Software
% C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
% C***********************************************************************

function [XX,YT,YC,XB,YB,NB,NAME] = naca4(IDES,NSIDE)
XX = zeros(1,NSIDE);
YT = zeros(1,NSIDE);
YC = zeros(1,NSIDE);
XB = zeros(1,2*NSIDE);
YB = zeros(1,2*NSIDE);

AN = 1.5;


N4 =  IDES                             / 1000;
N3 = (IDES - N4*1000                 ) / 100;
N2 = (IDES - N4*1000 - N3*100        ) / 10;
N1 = (IDES - N4*1000 - N3*100 - N2*10);

M = N4 / 100.0;
P = N3 / 10.0;
T = (N2*10 + N1) / 100.0;

ANP = AN + 1.0;
for ii = 1:NSIDE
    FRAC = (ii-1)/(NSIDE-1);
    if ii == NSIDE
        XX(ii) = 1.0;
    else
        XX(ii) = 1.0 - ANP*FRAC*(1.0-FRAC)^AN - (1.0-FRAC)^ANP;
    end
    YT(ii) = ( 0.29690*SQRT(XX(ii)) - 0.12600*XX(ii) - 0.35160*XX(ii)^2 + 0.28430*XX(ii)^3 - 0.10150*XX(ii)^4) * T / 0.20;
    if XX(ii) < P
        YC(ii) = M/P^2 * (2*P*XX(ii) - xx(ii)^2);
    else
        YC(ii) = M/(1-P)^2* ((1 -2*P) + 2*P*XX(ii) - XX(ii)^2);
    end
end

IB = 0;
for ii = NSIDE:-1:1
    IB = IB+1;
    XB(IB) = XX(ii);
    YB(IB) = YC(ii) + YT(ii);
end

for ii = 2:NSIDE
    IB = IB+1;
    XB(IB) = XX(ii);
    YB(IB) = YC(ii) - YT(ii);
end

NB = IB;
NAME = strcat('NACA',N4+1,N3+1,N2+1,N1+1);
