LIB "primdec.lib";
option(redSB);
ring r = 0,(a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p),dp;
poly p1 = pa-ek+fi-hm;
poly p2 = la-oa+bk+bp+ch-ci-dk-dp+el-eo-fg-fj+gm+hn-in+jm;
poly p3 = la+oa-bk-bp+ch+ci-dk-dp+el+eo-fg-fj-gm+hn+in-jm;
poly p4 = -ja+na+bf-bh-ce+cp+di-dm+eg-fo-gp+ho-il+jk-kn+lm;
poly p5 = ja+na+bf+bh-ce+cp-di-dm-eg-fo+gp-ho+il-jk-kn+lm;
poly p6 = -a-eka-epa+fia-jna+kpa-loa+bdk+bdp+bfg-bhn-cdi-ceg+cgp-cho+djm-e+ehm+ekp-elo-fip+fjo+glm-hkm-iln+jkn+k+p;
poly p7 = a-eka-epa+fia-jna+kpa-loa+bdk+bdp+bfg-bhn-cdi-ceg+cgp-cho+djm+e+ehm+ekp-elo-fip+fjo+glm-hkm-iln+jkn-k-p;
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

