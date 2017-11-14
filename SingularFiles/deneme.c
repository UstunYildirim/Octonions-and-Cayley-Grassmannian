LIB "primdec.lib";
option(redSB);
ring r = 0,(a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p),dp;
poly p1 = bl-dj-eo+gm;
poly p2 = ga-ja-bh+bi-ce-cl+df+dk+en-fm+gp-ho+io-jp-km+ln;
poly p3 = ga-ja+bh+bi-ce-cl-df+dk-en+fm-gp+ho+io+jp-km-ln;
poly p4 = gpa-hoa-jpa+lna-bgl+bhk+bip-blm-cep+cfl-chj+chm+d+deo-dfk+dgj-dgm-din+djm+ejo-ekn-fio+fkm-g+gin-gjm+j-m;
poly p5 = gpa-hoa-jpa+lna+bgl-bhk+bip-blm-cep-cfl+chj+chm+d+deo+dfk-dgj-dgm-din+djm-ejo+ekn+fio-fkm+g-gin+gjm-j-m;
poly p6 = fla-foa+gna-hja+b-bel+beo-bgm+bhi-bkp+blo-cen+cfm+cjp-cln+dej-dfi-djo+dkn-e+ekp-elo-gip+glm+hio-hkm-l+o;
poly p7 = -fla-foa+gna+hja+b+bel+beo-bgm-bhi-bkp+blo-cen+cfm+cjp-cln-dej+dfi-djo+dkn+e-ekp+elo+gip-glm-hio+hkm+l+o;
ideal I = p1,p2,p3,p4,p5,p6,p7;

I = groebner(I);
matrix H = jacob(I);
ideal J = I,(minor(H,4));

J = groebner(J);
J = radical(J);

print(J);
print(K);

