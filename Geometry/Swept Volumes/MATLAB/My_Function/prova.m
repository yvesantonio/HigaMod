clc
clear all
close all
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\model.stl');
stlPlot(vertices,faces,name);
open_verts = stlGetVerts(vertices,faces,'opened');
closed_verts = stlGetVerts(vertices,faces,'closed');
