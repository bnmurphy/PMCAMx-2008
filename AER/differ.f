       subroutine differ(rp,rl,fgl,flg,gfgl,gflg,dp,dl)
       dimension rp(28),rl(28),fgl(21),flg(21),gfgl(21),gflg(21)
       dimension dp(49),dl(49)
       do 1510 i=1,21
       dp(i)=rp(i)+fgl(i)
       dl(i)=rl(i)+flg(i)
1510   continue

       do 1520 i=22,28
       dp(i)=rp(i)
       dl(i)=rl(i)
1520   continue

       do 1530 i=29,49
       dp(i)=gflg(i-28)
       dl(i)=gfgl(i-28)
1530   continue


       return
       end
