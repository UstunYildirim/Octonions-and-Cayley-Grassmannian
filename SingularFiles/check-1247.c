LIB "primdec.lib";
option(redSB);
ring r = 0,(a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p),dp;
poly p1 = -bo+di+em-hk;
poly p2 = ma-oa+bl+bp+ch-ci-dg-dj-el-ep+fm-fo+gk+hn-in+jk;
poly p3 = ma-oa-bl+bp+ch+ci-dg-dj+el-ep-fm+fo-gk+hn+in-jk;
poly p4 = -fma+foa+ina-jka+b+bem-beo+bgn-bhk-blp+bmo-cej-cfh+cip-cjo+dei+dfg-d*im+djl-e+elp-emo-gkp+gmn+hko-hln+m-o;
poly p5 = fma-foa-ina+jka+b-bem+beo-bgn+bhk+blp-bmo+cej+cfh-cip+cjo-dei-dfg+d*im-djl-e-elp+emo+gkp-gmn-hko+hln+m-o;
poly p6 = ia-ka+bg-bn-ce-co-df+dl+ej+fh+gm-hl-ip+jo+kp-mn;
poly p7 = ia+ka+bg+bn+ce+co+df-dl+ej+fh+gm-hl-ip+jo-kp+mn;
ideal I = p1,p2,p3,p4,p5,p6,p7;

I = groebner(I);
matrix H = jacob(I);
ideal J = I,(minor(H,4));

J = groebner(J);
J = radical(J);


matrix M = jacob(J);
ideal K = J,(minor(M,11));
K = groebner (K);
print(J);
print(K);

