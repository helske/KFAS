t1<-ts(rnorm(100))
t2<-ts(rnorm(100))

t<-cbind(t1,t2)

model1<-SSModel(t~SSMcycle(period=10, type='common')+SSMcycle(period=10, type='distinct')+
                  SSMtrend(2,type="common")+SSMtrend(2,type="distinct")+
                  SSMseasonal(period=4,type="common")+SSMseasonal(period=4,type="distinct")+
                  SSMseasonal(period=4,type="common",sea.type="trig")+
                  SSMseasonal(period=4,type="distinct",sea.type="trig")
                )
