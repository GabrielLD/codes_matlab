function [checkerBoard] = checkerboardMR(size_x,size_y,square_Side)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%sx=1280;sy=960;
%bas=16     % length base square side
% assume sx and sy are multiples of base

Lx=(-1).^[1:1:size_x/square_Side];
Ly=(-1).^[1:1:size_y/square_Side];

A=ones(size_x,size_y); %generate an array of zeros the size of the final image
[Ai Aj]=ind2sub([size_x size_y],[1:1:size_x*size_y]); % construction of the indices of the elements in matrix A
Ai2=reshape(Ai,[size_x size_y]);
Aj2=reshape(Aj,[size_x size_y]);
linex=[1:square_Side:size_x]; %position of the top side of squares
linex=[linex size_x];
liney=[1:square_Side:size_y];
liney=[liney size_y];

setx=zeros(1,square_Side);
for k=2:1:(length(linex)-1)
  L1=[linex(k-1):1:(linex(k)-1)]
       setx=[setx;L1];  
end
setx(1,:)=[];
sety=zeros(1,square_Side);
for k=2:1:length(liney)-1 
  L2=[liney(k-1):1:liney(k)-1];
  sety=[sety;L2];  
end
sety(1,:)=[];
blk=ones(square_Side,square_Side);
for k=1:1:size_x/square_Side-1
  for s=1:1:size_y/square_Side-1
          Lx1=blk*Lx(k);
            Ly1=blk*Ly(s);
            L12=Lx1.*Ly1;
          A(setx(k,:),sety(s,:))=A(setx(k,:),sety(s,:)).*L12;
%             A(sety(s,:),setx(k,:))=A(sety(s,:),setx(k,:)).*L12;
  end;  
end;
;
checkerBoard=A;
end

