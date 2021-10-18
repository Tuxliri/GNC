function dX = f(t,X,mu)

% Function providing dx(t)=f(x,t)
x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);

r2 = x^2+y^2+z^2;
dx = [vx;
      vy;
      vz;
      -mu*x/(r2)^(3/2);
      -mu*y/(r2)^(3/2);
      -mu*z/(r2)^(3/2)];
  
  %% MAYBE THERE IS SOMETHING WRONG WITH THE STM DEFINITION!
  PHI11=X(7);
  PHI12=X(8);
  PHI13=X(9);
  PHI14=X(10);
  PHI15=X(11);
  PHI16=X(12);
  PHI21=X(13);
  PHI22=X(14);
  PHI23=X(15);
  PHI24=X(16);
  PHI25=X(17);
  PHI26=X(18);
  PHI31=X(19);
  PHI32=X(20);
  PHI33=X(21);
  PHI34=X(22);
  PHI35=X(23);
  PHI36=X(24);
  PHI41=X(25);
  PHI42=X(26);
  PHI43=X(27);
  PHI44=X(28);
  PHI45=X(29);
  PHI46=X(30);
  PHI51=X(31);
  PHI52=X(32);
  PHI53=X(33);
  PHI54=X(34);
  PHI55=X(35);
  PHI56=X(36);
  PHI61=X(37);
  PHI62=X(38);
  PHI63=X(39);
  PHI64=X(40);
  PHI65=X(41);
  PHI66=X(42);
  
dPHI=[PHI41;
      PHI42;
      PHI43;
      PHI44;
      PHI45;
      PHI46;
      PHI51;
      PHI52;
      PHI53;
      PHI54;
      PHI55;
      PHI56;
      PHI61;
      PHI62;
      PHI63;
      PHI64;
      PHI65;
      PHI66;
      (3*PHI21*mu*x*y)/(r2)^(5/2) - PHI11*(mu/(r2)^(3/2) - (3*mu*x^2)/(r2)^(5/2)) + (3*PHI31*mu*x*z)/(r2)^(5/2);
      (3*PHI22*mu*x*y)/(r2)^(5/2) - PHI12*(mu/(r2)^(3/2) - (3*mu*x^2)/(r2)^(5/2)) + (3*PHI32*mu*x*z)/(r2)^(5/2);
      (3*PHI23*mu*x*y)/(r2)^(5/2) - PHI13*(mu/(r2)^(3/2) - (3*mu*x^2)/(r2)^(5/2)) + (3*PHI33*mu*x*z)/(r2)^(5/2);
      (3*PHI24*mu*x*y)/(r2)^(5/2) - PHI14*(mu/(r2)^(3/2) - (3*mu*x^2)/(r2)^(5/2)) + (3*PHI34*mu*x*z)/(r2)^(5/2);
      (3*PHI25*mu*x*y)/(r2)^(5/2) - PHI15*(mu/(r2)^(3/2) - (3*mu*x^2)/(r2)^(5/2)) + (3*PHI35*mu*x*z)/(r2)^(5/2);
      (3*PHI26*mu*x*y)/(r2)^(5/2) - PHI16*(mu/(r2)^(3/2) - (3*mu*x^2)/(r2)^(5/2)) + (3*PHI36*mu*x*z)/(r2)^(5/2);
      (3*PHI11*mu*x*y)/(r2)^(5/2) - PHI21*(mu/(r2)^(3/2) - (3*mu*y^2)/(r2)^(5/2)) + (3*PHI31*mu*y*z)/(r2)^(5/2);
      (3*PHI12*mu*x*y)/(r2)^(5/2) - PHI22*(mu/(r2)^(3/2) - (3*mu*y^2)/(r2)^(5/2)) + (3*PHI32*mu*y*z)/(r2)^(5/2);
      (3*PHI13*mu*x*y)/(r2)^(5/2) - PHI23*(mu/(r2)^(3/2) - (3*mu*y^2)/(r2)^(5/2)) + (3*PHI33*mu*y*z)/(r2)^(5/2);
      (3*PHI14*mu*x*y)/(r2)^(5/2) - PHI24*(mu/(r2)^(3/2) - (3*mu*y^2)/(r2)^(5/2)) + (3*PHI34*mu*y*z)/(r2)^(5/2);
      (3*PHI15*mu*x*y)/(r2)^(5/2) - PHI25*(mu/(r2)^(3/2) - (3*mu*y^2)/(r2)^(5/2)) + (3*PHI35*mu*y*z)/(r2)^(5/2);
      (3*PHI16*mu*x*y)/(r2)^(5/2) - PHI26*(mu/(r2)^(3/2) - (3*mu*y^2)/(r2)^(5/2)) + (3*PHI36*mu*y*z)/(r2)^(5/2);
      (3*PHI11*mu*x*z)/(r2)^(5/2) - PHI31*(mu/(r2)^(3/2) - (3*mu*z^2)/(r2)^(5/2)) + (3*PHI21*mu*y*z)/(r2)^(5/2);
      (3*PHI12*mu*x*z)/(r2)^(5/2) - PHI32*(mu/(r2)^(3/2) - (3*mu*z^2)/(r2)^(5/2)) + (3*PHI22*mu*y*z)/(r2)^(5/2);
      (3*PHI13*mu*x*z)/(r2)^(5/2) - PHI33*(mu/(r2)^(3/2) - (3*mu*z^2)/(r2)^(5/2)) + (3*PHI23*mu*y*z)/(r2)^(5/2);
      (3*PHI14*mu*x*z)/(r2)^(5/2) - PHI34*(mu/(r2)^(3/2) - (3*mu*z^2)/(r2)^(5/2)) + (3*PHI24*mu*y*z)/(r2)^(5/2);
      (3*PHI15*mu*x*z)/(r2)^(5/2) - PHI35*(mu/(r2)^(3/2) - (3*mu*z^2)/(r2)^(5/2)) + (3*PHI25*mu*y*z)/(r2)^(5/2);
      (3*PHI16*mu*x*z)/(r2)^(5/2) - PHI36*(mu/(r2)^(3/2) - (3*mu*z^2)/(r2)^(5/2)) + (3*PHI26*mu*y*z)/(r2)^(5/2)];
  
dX = [dx;dPHI];
end

