% %
%     Copyright (C) 2024, Andria Pelentritou
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% %
%

function [ arear ] = rocarea2(ff,tt)

ttn = [1, tt];
ffn = [1, ff];
f_step = abs(diff(ffn));
hd1=ttn(1,1:end-1).*f_step;
arear1 = sum(hd1);
hd2=ttn(1,2:end).*f_step;
arear2 = sum(hd2);
arear=(arear1+arear2)/2;