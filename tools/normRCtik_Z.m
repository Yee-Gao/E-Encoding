 function [h,CR_test]=normRCtik_Z(x,y,D,tikf)
%normRCtik_Z    ����TRF
%�ɽ���������������x��y�ֱ�Ϊ���롢�������ÿ�ж�Ӧһ������������
%D��TRF������tikfΪ���򻯲�����hΪTRFģ�ͣ�CR_tesΪģ�͵�Ԥ��������
% x, stimulus matrix; y, response matrix; each stimulus/response is a row vector
% D, the order of the TRF;
% tikf, regularization parameter;
h=zeros(D,size(x,1),size(y,1),10);
CR_test=zeros(size(y,1),10);
for segno=0:9               %��ѭ������ʵ�ֽ�����飬ÿ��ѭ����ȡһ����������Ϊ��������
    [h(:,:,:,segno+1),CR_test(:,segno+1)]=normRCtik(x,y,D,segno,tikf);
end
%h����=����������=�����źŸ���������ά��=����źŸ���
h=mean(h,4);
CR_test=mean(CR_test,2);

function [h,CR_test,Str_testE,CR_train]=normRCtik(x,y,D,segno,tikf)
% NRC with Tikhonov regularization
% each stimulus/response is a row vector
CR_test=zeros(size(y,1),1);
TSlen=size(x,2)/10; %testing sample number
TSlen=floor(TSlen);
%vector form
testing_range=[1:TSlen]+TSlen*segno;
training_range=setdiff([1:length(x)],testing_range);
x_test=x(:,testing_range);
y_test=y(:,testing_range);
x=x(:,training_range);
y=y(:,training_range);
dimx=size(x,1);
lenx=size(x,2);
order=D;
chhn=size(y,1);
Y=y(:,order:lenx);           %����ȥ��ǰ�ˣ�D-1����Ԫ��
%% 
X=zeros(dimx*order,lenx-order+1);
for ind1=1:dimx
  for ind2=1:order
    X((ind1-1)*order+ind2,:)=x(ind1,ind2:(lenx-order+ind2));
  end
end
%% 
if tikf==Inf
  h=X*Y';
else
  Rx=X*X';
  d=eig(Rx);
  tikf=tikf*d(end);
  h=(Rx+tikf*eye(size(Rx)))\(X*Y');
end
h=reshape(h,order,dimx,chhn);
h=h(end:-1:1,:,:);              %�˴�h����Ϊ����������Ϊ�����źŸ���
%%
y_pred=y_test*0;
for chh=1:chhn
    for ind1=1:size(h,2)
      y_pred(chh,:)=y_pred(chh,:)+filter(h(:,ind1,chh),1,x_test(ind1,:));
    end
    CR_test(chh,:)=sum(y_test(chh,:).*y_pred(chh,:))/sqrt(sum(y_pred(chh,:).*y_pred(chh,:))*sum(y_test(chh,:).*y_test(chh,:)));
end
