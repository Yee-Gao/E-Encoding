 function [h,CR_test]=normRCtik_Z(x,y,D,tikf)
%normRCtik_Z    计算TRF
%可解多输入多输出情况；x、y分别为输入、输出矩阵，每行对应一个输入或输出；
%D是TRF阶数，tikf为正则化参数，h为TRF模型，CR_tes为模型的预测能力；
% x, stimulus matrix; y, response matrix; each stimulus/response is a row vector
% D, the order of the TRF;
% tikf, regularization parameter;
h=zeros(D,size(x,1),size(y,1),10);
CR_test=zeros(size(y,1),10);
for segno=0:9               %此循环用来实现交叉检验，每次循环抽取一部分数据做为测试数据
    [h(:,:,:,segno+1),CR_test(:,segno+1)]=normRCtik(x,y,D,segno,tikf);
end
%h行数=阶数，列数=输入信号个数，第三维数=输出信号个数
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
Y=y(:,order:lenx);           %各行去除前端（D-1）个元素
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
h=h(end:-1:1,:,:);              %此处h行数为阶数，列数为输入信号个数
%%
y_pred=y_test*0;
for chh=1:chhn
    for ind1=1:size(h,2)
      y_pred(chh,:)=y_pred(chh,:)+filter(h(:,ind1,chh),1,x_test(ind1,:));
    end
    CR_test(chh,:)=sum(y_test(chh,:).*y_pred(chh,:))/sqrt(sum(y_pred(chh,:).*y_pred(chh,:))*sum(y_test(chh,:).*y_test(chh,:)));
end
