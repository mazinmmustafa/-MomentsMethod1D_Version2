close all; clear; clc;

R = 1;

Basis = load("basis.dat");

[N,~] = size(Basis);

figure()

##view(45, 30)
view(0, 0)

axis equal
##axis([-1 +1 -1 +1 -1 +1]*R)
##xlim([-1 +1])
##xlim([+1.1 +1.3])
##zlim([-0.5 +0.5])

hold on

for i=1:N
	if Basis(i,10)==0
		plot3([Basis(i,1), Basis(i,4)], [Basis(i,2), Basis(i,5)], [Basis(i,3), Basis(i,6)],'b','LineWidth',3)
		plot3([Basis(i,4), Basis(i,7)], [Basis(i,5), Basis(i,8)], [Basis(i,6), Basis(i,9)],'r','LineWidth',3)
	else
		plot3([Basis(i,1), Basis(i,4)], [Basis(i,2), Basis(i,5)], [Basis(i,3), Basis(i,6)],'-g','LineWidth',8)
		plot3([Basis(i,4), Basis(i,7)], [Basis(i,5), Basis(i,8)], [Basis(i,6), Basis(i,9)],'-c','LineWidth',8)
	end
	input("Next");
end

hold off




