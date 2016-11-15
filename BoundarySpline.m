% This scripts recieves the initial point to point path from a planner or
% user. Several algorithms are run (Midpoint, Sharp, Heading I, Heading F)
% The resulting path is plotted using a non-uniform clamped B-spline.
% B-spline curve is evaluated using the Cox-deBoor recursvie algorithm.
%                           
%
%
% Written by Mohamed Elbanhawi
% Date: 06/11/14
%
% Completed: Midpoints, Sharpness, Initial and Final Headings B-Spline Plot
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BoundarySpline(x0,y0,t0,xf,yf,tf)
    %Setup figure
    figure (1);
    
    %Segment Parameters
    global  Lmin;
    global Amin;
    Lmin = 7.5;
    Amin =pi/2;
    
    %Control Polylines
    px = [];
    py = [];
  
    x1 = x0+Lmin*cos(t0);
    y1 = y0+Lmin*sin(t0);
    
    x2 = xf+Lmin*cos(tf-pi);
    y2 = yf+Lmin*sin(tf-pi); 
    
    px=[x0,x1,x2,xf];
    py=[y0,y1,y2,yf];

    %Ensure that all turns do not exceed turning angle thershold
    
    %--------Start segment--------%
    %check start
    [sx,sy]=sharpb(px,py);
    
    
    px=[sx,px(end-1:end)];
    py=[sy,py(end-1:end)];

    %--------End Segment--------%
    [ex,ey]=sharpb(fliplr(px),fliplr(py));
    ex =  fliplr(ex);
    ey =  fliplr(ey);
    
    %%update control polyline
    px = [sx,ex];
    py = [sy,ey];
     
    
    %Add Midpoints
    [px,py]=midpoint(px,py);
   
    
    %check mid section feasbility
    [px,py]=sharp(px,py);
    
    %Produce plots
    [x,y]=Bspline(px,py);
    
    %Curvature(x,y,0.1);
    %Tracking(x,y);
    %HeadingTrack(x,y,0.1,px(1),py(1),qtan(py(1),px(1),py(2),px(2)),10000);

  
end

function [px,py]=midpoint(px,py)
    %add midpoints
    for i=0:1:size(px,2)-1;
        mx(2*i+1)=px(i+1);
        my(2*i+1)=py(i+1);
    end
    for i=2:2:size(mx,2)-1;
        mx(i)=(mx(i-1)+mx(i+1))/2;
        my(i)=(my(i-1)+my(i+1))/2;
    end
    px=mx(:)';
    py=my(:)';
end

function [px,py]=sharp(px,py)
    global  Lmin;
    global  Amin;
    %number of segments = (number of points-5)/2 + 1
    %seg=(size(px,2)-3)/2;
    i=2;
    while i <((size(px,2)-3)/2);
        n=2*i-1; %segment start point
        
        %Apply cosine rule to find turning angle value
        c=sqrt((px(n)-px(n+4))^2+(py(n)-py(n+4))^2);
        a=sqrt((px(n)-px(n+2))^2+(py(n)-py(n+2))^2);
        b=sqrt((px(n+2)-px(n+4))^2+(py(n+2)-py(n+4))^2);
        H=(a*a+b*b-c*c)/(2*a*b);
        gamma= acos(H);
        gamma*180/pi;    
       
        
        if gamma<(Amin)
            %shift cells in arrays by two starting from the last cell in
            %segment
            for j=size(px,2):-1:n+4;
                px(j+2)=px(j);
                py(j+2)=py(j);
            end
%             %calculate the three new points
%             %locate points relative to each other
%             s(1,1)=px(n);s(1,2)=py(n);
%             g(1,1)=px(n+2);g(1,2)=py(n+2);
%             v(1,1)=px(n+6);v(1,2)=py(n+6);
%             d= sign(det([(s-g);(g-v)])/(norm(s-g)));
%             dx=sign(px(n+6)-px(n))*1;
%             dy=sign(py(n+6)-py(n))*1;
            
            %New points
            ang=qtan(py(n),px(n),py(n+2),px(n+2));
            
            px1=px(n+2)+Lmin*cos(pi/2+ang);
            py1=py(n+2)+Lmin*sin(pi/2+ang);
            d1 = sqrt((px1-px(n+6))^2+(py1-py(n+6))^2);
            px2=px(n+2)-Lmin*cos(pi/2+ang);
            py2=py(n+2)-Lmin*sin(pi/2+ang);
            d2 = sqrt((px2-px(n+6))^2+(py2-py(n+6))^2);
            
            if d1<=d2;
                px(n+4) = px1;
                py(n+4) = py1;
            else
                px(n+4) = px2;
                py(n+4) = py2;
            end
            
            px(n+5)=(px(n+4)+px(n+6))/2;
            py(n+5)=(py(n+4)+py(n+6))/2;
            
            px(n+3)=(px(n+4)+px(n+2))/2;
            py(n+3)=(py(n+4)+py(n+2))/2;
            i=i+1;
        end
       
       i=i+1;
       
    end
end

function [px,py]=sharpb(px,py)
    global  Lmin;
    global  Amin;
   
    i=1;
    
        n=i; %segment start point
        
        %Apply cosine rule to find turning angle value
        c=sqrt((px(n)-px(n+2))^2+(py(n)-py(n+2))^2);
        a=sqrt((px(n)-px(n+1))^2+(py(n)-py(n+1))^2);
        b=sqrt((px(n+1)-px(n+2))^2+(py(n+1)-py(n+2))^2);
        H=(a*a+b*b-c*c)/(2*a*b);
        gamma= acos(H);
        
        if gamma<(Amin)
            %shift cells in arrays by two starting from the last cell in
            %segment
            for j=size(px,2):-1:n+2;
                px(j+1)=px(j);
                py(j+1)=py(j);
            end
            %calculate the three new points
            %locate points relative to each other
%             s(1,1)=px(n);s(1,2)=py(n);
%             g(1,1)=px(n+1);g(1,2)=py(n+1);
%             v(1,1)=px(n+2);v(1,2)=py(n+2);
%             d= sign(det([(s-g);(g-v)])/(norm(s-g)));
%             
%      
%             dx=sign(px(n+3)-px(n))*1;
%             dy=sign(py(n+3)-py(n))*1;
%            
            %New points
            ang=qtan(py(n),px(n),py(n+1),px(n+1));
            
            px1=px(n+1)+Lmin*cos(pi/2+ang);
            py1=py(n+1)+Lmin*sin(pi/2+ang);
            d1 = sqrt((px1-px(n+3))^2+(py1-py(n+3))^2);
            px2=px(n+1)-Lmin*cos(pi/2+ang);
            py2=py(n+1)-Lmin*sin(pi/2+ang);
            d2 = sqrt((px2-px(n+3))^2+(py2-py(n+3))^2);
            
            if d1<=d2;
                px(n+2) = px1;
                py(n+2) = py1;
            else
                px(n+2) = px2;
                py(n+2) = py2;
            end
            
            px =px(1:3);
            py =py(1:3);
        else
            px =px(1:2);
            py =py(1:2);
        end
       
       i=i+1;
       
    end

function [x,y]=Bspline(px,py)

px=px(:)';py=py(:)';
%Intialize some variables
du = 0.001;
u=0:du:1; %parameter
x=u*0;
y=u*0;
% Calculate basis function for a B-spline curve based on user inputed control points 
for deg=3;
% generate normalized uniform knot vector with clamped ends
uhat=knot(px,deg);
% evaluate deBoor algorithm at each point along the parameter values
for i=1:1:size(u,2);
    [x(i),y(i)]=deBoor(u(i),px,py,uhat,deg);
end
figure(1); plot(x,y,'LineWidth',2, 'Color',[34,147,216]/255);grid off; xlabel ('x','FontName','Times New Roman','FontAngle','Italic','FontSize',22); ylabel('y','FontName','Times New Roman','FontAngle','Italic','FontSize',22);hold on;
%figure(1); plot(px,py,'-*','Color',[165,165,165]/255);hold on
axis('equal')
end
end

function [x,y]=deBoor(u,px,py,uhat,deg)
%create basis function empty matrix
N=zeros(size(px,2)+deg,deg+1); % degree+1 = 4th order is assumed for cubic curves
%create first level of basis function where p=1, and i<= control+3
for i=1:1:size(px,2)+deg;
    if ((u >= uhat(i)) && (u < uhat(i+1)));
        N(i,1)=1;
    end
end
if u==1
    for i=1:1:size(px,2)+deg;
    if ((u >= uhat(i)) && (u <= uhat(i+1)));
        N(i,1)=1;
    end
    end
end
%de Boor's recursive algorithm
for p=2:1:deg+1;
    i_max = size(px,2)+deg+1-p;
    for i=1:1:i_max;
        if uhat(i+p-1)~=uhat(i) && uhat(i+p)~=uhat(i+1)
        N(i,p)=(N(i,p-1)*(u-uhat(i))/(uhat(i+p-1)-uhat(i)))+((N(i+1,p-1)*(-u+uhat(i+p))/(uhat(i+p)-uhat(i+1))));
        %N(isnan(N)) = 0 ;
        elseif uhat(i+p-1)==uhat(i)
        N(i,p)=((N(i+1,p-1)*(-u+uhat(i+p))/(uhat(i+p)-uhat(i+1))));
        elseif uhat(i+p)==uhat(i+1)
        N(i,p)=(N(i,p-1)*(u-uhat(i))/(uhat(i+p-1)-uhat(i)));
        end
    end
end
%Evaluate bspline curve at this instance
x=0;
y=0;
for i=1:1:size(px,2)
    x=x+N(i,deg+1)*px(i);
    y=y+N(i,deg+1)*py(i);
end
end

function uhat=knot(px,p)
knot_size=size(px,2)+(p+1); % m=n+p+1
delta_uhat=1/(knot_size+1-(p+1)*2); %end knots are of the orders multiplicity
uhat=zeros(1,knot_size);
for i=1:1:p+1
    uhat(i)=0;
    uhat(knot_size+1-i)=1;
end

for i=p+2:knot_size-(p+1)
    uhat(i)=uhat(i-1)+delta_uhat;
end
end

function theta=qtan(y1,x1,y2,x2)
    
theta=atan((y2-y1)/(x2-x1));

%First Quadrant
if (x2>x1)&&(y2>y1);
    theta =abs(theta);
    
%90 degrees
elseif (x2==x1)&&(y2>y1);
    theta = pi/2;
   
%180 degrees
elseif (x2<x1)&&(y2==y1);
    theta = pi;
    
%270 degrees
elseif (x2==x1)&&(y2<y1);
    theta =-pi/2;   
    
%Second Quadrant
elseif (x2<x1)&&(y2>y1);
    theta = pi - abs(theta);
    
%Third Quadrant
elseif (x2<x1)&&(y2<y1);
    theta = abs(theta) -pi;
    
%Fourth Quadrant
elseif (x2>x1)&&(y2<y1);
    theta = -abs(theta);
    
%Zero
elseif (x2>x1)&&(y2==y1);
    theta =0;
    
end

end


