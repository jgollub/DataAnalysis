function draw_extent_box(Azex,Elex,Zex,obj_numbers,color,alpha)

if nargin==6
    a = alpha;
else
    a = 0.1;
end

if nargin==5
    c = color;
else
    c = [0,1,0.2];
end

for n=obj_numbers
    
    x1 = tan(Azex(n,1))*Zex(n,1);
    x2 = x1;
    x3 = tan(Azex(n,2))*Zex(n,1);
    x4 = x3;
    x5 = tan(Azex(n,1))*Zex(n,2);
    x6 = x5;
    x7 = tan(Azex(n,2))*Zex(n,2);
    x8 = x7;

    y1 = tan(Elex(n,1))*Zex(n,1)*cos(Azex(n,1));
    y2 = tan(Elex(n,2))*Zex(n,1)*cos(Azex(n,1));
    y3 = y2;
    y4 = y1;
    y5 = tan(Elex(n,1))*Zex(n,2)*cos(Azex(n,1));
    y6 = tan(Elex(n,2))*Zex(n,2)*cos(Azex(n,1));
    y7 = y6;
    y8 = y5;

    z1 = Zex(n,1);
    z2 = Zex(n,1);
    z3 = Zex(n,1);
    z4 = Zex(n,1);
    z5 = Zex(n,2);
    z6 = Zex(n,2);
    z7 = Zex(n,2);
    z8 = Zex(n,2);

    patch([x1 x2 x3 x4], [z1 z2 z3 z4], [y1 y2 y3 y4],c,'FaceAlpha',a)
    patch([x5 x6 x7 x8], [z5 z6 z7 z8], [y5 y6 y7 y8],c,'FaceAlpha',a)
    patch([x1 x5 x8 x4], [z1 z5 z8 z4], [y1 y5 y8 y4],c,'FaceAlpha',a)
    patch([x2 x6 x7 x3], [z2 z6 z7 z3], [y2 y6 y7 y3],c,'FaceAlpha',a)
    patch([x2 x6 x5 x1], [z2 z6 z5 z1], [y2 y6 y5 y1],c,'FaceAlpha',a)
    patch([x2 x6 x5 x1], [z2 z6 z5 z1], [y2 y6 y5 y1],c,'FaceAlpha',a)
    patch([x3 x7 x8 x4], [z3 z7 z8 z4], [y3 y7 y8 y4],c,'FaceAlpha',a)
end