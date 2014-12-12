function proceding_region=findWhichLPlotQuarter(data1,control_var1,data2,control_var2)

data_bin(:,1)=data1<control_var1;
data_bin(:,2)=data2<control_var2;

data_xor=xor(data_bin(:,1),data_bin(:,2));
indx=find(data_xor==1);
diff=data_bin(indx,1)-data_bin(indx,2);
nr_region1=length(find(diff==-1))
nr_region2=length(find(diff==1))

if nr_region1>nr_region2
    proceding_region='region1';
else
    proceding_region='region2';
end