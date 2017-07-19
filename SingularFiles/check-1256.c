LIB "primdec.lib";
option(redSB);
ring r = 0,(a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p),dp;
poly p1 = -pa+di+fl-hk;
poly p2 = ma-oa+bl+bp+ch-ci-dg-dj-el-ep+fm-fo+gk+hn-in+jk;
poly p3 = -ma-oa+bl+bp+ch+ci-dg-dj+el+ep-fm-fo-gk+hn+in-jk;
poly p4 = -ja+na+bi-bk-cf+cp+de-dm-eh+fg-gp+hm-io+jl+ko-ln;
poly p5 = -ja-na+bi+bk+cf-cp-de+dm-eh+fg-gp+hm-io+jl-ko+ln;
poly p6 = a+fla+fpa-hka+jna-lpa+moa-bel-bep+bgk-bin-ceh+cfg-cgp-cio+dej-dfi+dgm+dil+f-flp+fmo+hkp-hmn+jko-jln-l-p;
poly p7 = a-fla-fpa+hka-jna+lpa-moa+bel+bep-bgk+bin+ceh-cfg+cgp+cio-dej+dfi-dgm-dil+f+flp-fmo-hkp+hmn-jko+jln-l-p;
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

