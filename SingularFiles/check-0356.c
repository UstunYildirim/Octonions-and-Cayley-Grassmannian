LIB "primdec.lib";
option(redSB);
ring r = 0,(a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p),dp;
poly p1 = bo-dl+fi-hm;
poly p2 = la-oa+bk+bp+ch-ci-dk-dp+el-eo-fg-fj+gm+hn-in+jm;
poly p3 = -la+oa+bk-bp+ch+ci-dk+dp+el-eo-fg-fj-gm+hn+in-jm;
poly p4 = ela-eoa-fja+hna-b-bdl+bdo+bfi-bgn+bkp-blo+cdj-cei+chp+cjo+d-dhm-dkp+dlo+egm-fgp-fio-gln+hlm-ikn+jkm-l+o;
poly p5 = ela-eoa-fja+hna+b-bdl+bdo+bfi-bgn+bkp-blo+cdj-cei+chp+cjo-d-dhm-dkp+dlo+egm-fgp-fio-gln+hlm-ikn+jkm+l-o;
poly p6 = fa-ha+bg-bn-cd-co+dj-ei+em-fp+gl+hp+ik+jo-km-ln;
poly p7 = fa+ha-bg-bn-cd-co-dj+ei+em-fp-gl-hp-ik-jo-km-ln;
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

