Name = 'T0803'
leftfile = ['/home/xun/Downloads/',Name,'-contact.txt']
rightfile = ['/home/xun/Downloads/',Name,'-contactmap.txt']
[A,delimit1] = importdata(leftfile);
[B,delimit2] = importdata(rightfile);
x = A(:,3);
y = B(:,3);
for i = 1 : length(x);
   if A(i,3) == 1
    rectangle('Position',[A(i,1),A(i,2),1,1],'Facecolor',[0,0,A(i,3)])
   end
end
for i = 1 : length(y)
   if B(i,3) > 0.5
     rectangle('Position',[B(i,2),B(i,1),1,1],'Facecolor',[1,0,0])
   end
   %if B(i,3) < 0.5
    % rectangle('Position',[B(i,2),B(i,1),1,1],'Facecolor',[1,1,0])
   %end
end
ylabel('native')
xlabel('contactmap')
title([Name,' contactmap blue - native red contact > 0.5'])
saveas(gcf,['/home/xun/data_analysis/',Name,'-contactmap0.5.png'])

