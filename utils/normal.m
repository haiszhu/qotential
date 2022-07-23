function Nc = normal(A,B,C)
% normal to countclockwise orientation triangular surface
Nc=cross(B-A,C-A);
Nc = Nc(:)/norm(Nc);