 function im_pad_temp = pcalculation(iii,Dr,elmspace,PAD)
     im_pad_temp=zeros(2*Dr+1,2*Dr+1,2*Dr+1);
      for ii=3:30
         
    xi=Dr+1+(-16.5+iii)*elmspace;
    yi=Dr+1+(-16.5+ii)*elmspace;
    zi=1;
    for x=1:2*Dr+1
        for y=1:2*Dr+1
            for z=1:2*Dr+1
%             if sqrt((x-Dr-1)^2+(y-Dr-1)^2)<Dr
                r=floor((sqrt((x-xi)^2+(y-yi)^2+(z-zi)^2)));
                ppi=r-0;
                if ppi<=size(PAD,1)-1 && ppi>20 && z/r>0.997
%                     im_pad(x,y)=im_pad(x,y)+(PA(ppi,iii)-r*(PA(ppi+1,iii)-PA(ppi,iii)))/y/tl(r);
%                     im_pad(x,y)=im_pad(x,y)+(PA(ppi,iii))*y;
%                         im_pad(x,y)=im_pad(x,y)+(PA(ppi,iii))*y/r.^1.2;
%                         im_pad2(x,y)=im_pad2(x,y)+(PA2(ppi,iii))*y/r.^1.2;
                        im_pad_temp(x,y,z)=im_pad_temp(x,y,z)+PAD(ppi,iii,ii)*z/r.^2;
%                     im_pad_dump(x,y)=im_pad_dump(x,y)+1;
%                     im_pad_temp(x,y,z)=im_pad_temp(x,y,z)+(PAD(ppi,iii,ii)-(1.54./(0.7812*4))^2*r*(PAD(ppi+1,iii,ii)-PAD(ppi,iii,ii)));
%                    im_pad_temp(x,y,z)=im_pad_temp(x,y,z)+PAD(ppi,iii,ii);
                end
             end
            
        end
%         disp(x);
    end
      end