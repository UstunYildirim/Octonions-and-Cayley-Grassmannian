LIB "primdec.lib";
option(redSB);
ring r = 0,(a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p),dp;
poly p1 = -cn+fi+gh-jm;
poly p2 = hla-ika-joa+kna+bfl+bgk+bjp-bln+c-cdo+cep-cfi+cfn-cgh+cin-dhp+dio-dkm-efp-ego+elm-f-fin+fjm-ghn-i+ijm+n;
poly p3 = -hla+ika+joa-kna-bfl-bgk-bjp+bln+c+cdo-cep+cfi-cfn+cgh-cin+dhp-dio+dkm+efp+ego-elm-f+fin-fjm+ghn-i-ijm+n;
poly p4 = -ia+na+bg-bj-ce+cp+dh-dm+ef-fp-go-hl+ik+jo-kn+lm;
poly p5 = -ia+na+bg+bj+ce+cp-dh-dm-ef-fp-go+hl-ik-jo+kn+lm;
poly p6 = ha-ja+bf-bn-cd-co+di-eg+em+fl+gk-hp+io+jp-km-ln;
poly p7 = ha+ja+bf-bn+cd-co-di+eg+em-fl-gk-hp+io-jp-km+ln;
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

