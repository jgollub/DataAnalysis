y_length=66;
x_length=104;

upsample=4;
NSI_Prop_fields_complex=[];
for i=1:y_length
NSI_Prop_fields_abs(i,1:x_length)=NSI_Prop_fields(x_length*(i-1)+1:i*x_length,3);
NSI_Prop_fields_phase(i,1:x_length)=NSI_Prop_fields(x_length*(i-1)+1:i*x_length,4);
NSI_Cordinates_X(i,1:x_length)=NSI_Prop_fields(x_length*(i-1)+1:i*x_length,1);
NSI_Cordinates_Y(i,1:x_length)=NSI_Prop_fields(x_length*(i-1)+1:i*x_length,2);
end
NSI_RE_IM=10.^(NSI_Prop_fields_abs./20).*exp(1i*NSI_Prop_fields_phase*pi/180);
%imagesc(NSI_Cordinates_X(1,:),NSI_Cordinates_Y(:,1),flipdim(upsample_image(NSI_Prop_fields_abs,1),1),[-30,0]);

figure(1)
norm_NSI=max(max(NSI_RE_IM))
%imagesc(NSI_Cordinates_X(1,:),NSI_Cordinates_Y(:,1),flipdim(upsample_image(NSI_Prop_fields_abs,4),1),[-30,0]);
%imagesc(NSI_Cordinates_X(1,:),NSI_Cordinates_Y(:,1),flipdim(upsample_image(NSI_Prop_fields_phase,1),1));
 
imagesc(NSI_Cordinates_X(1,:),NSI_Cordinates_Y(:,1),upsample_image(abs(fliplr(flipud(NSI_RE_IM./norm_NSI))),upsample));
% imagesc(NSI_Cordinates_X(1,:),NSI_Cordinates_Y(:,1),upsample_image(angle(fliplr(flipud(NSI_RE_IM))),upsample));


axis equal
axis tight
 axis([-.25,.25,-.25,.25])
colorbar
% load('C:\Users\Jonah Gollub\Downloads\A4_Pol_M_NFS-23-Nov-2013 (1).mat')
% Propagate_Single_Field_Slice
figure(2)
norm=max(max(abs(Ey_Pan(:,:,1,51,6))))
imagesc(Az(1,:),El(:,1),flipdim(upsample_image(abs(Ey_Pan(:,:,1,51,6)./norm),upsample),2))
% imagesc(Az(1,:),El(:,1),flipdim(upsample_image(angle(Ey_Pan(:,:,1,51,6)),upsample),2))
axis equal
axis tight
axis([-.25,.25,-.25,.25])
colorbar

figure(3)

NSI_Subtract=(NSI_RE_IM./norm_NSI-Ey_Pan(:,:,1,51,6)./norm-NSI_RE_IM/norm_NSI);
 imagesc(NSI_Cordinates_X(1,:),NSI_Cordinates_Y(:,1),upsample_image(abs(NSI_Subtract),upsample))
% imagesc(NSI_Cordinates_X(1,:),NSI_Cordinates_Y(:,1),upsample_image(angle(NSI_Subtract),upsample))
axis equal
axis tight
colorbar
% figure(4)
%plot(Ey_Pan(33,:,1,51,6)
