function fXis_1(nb,TiempoC,BM)
      Ints=BM
      Xis=Ints
      return Xis
end
##
function fXis_2(nb,TiempoC,BM)
      Ints=BM
       Xis=sqrt(2)*H2(Ints)
      return Xis 
end
##
function fXis_3(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(3))*H3(Ints)
      return Xis
end
##
function fXis_4(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(4))*H4(Ints)
      return Xis
end
#   
function fXis_5(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(5))*H5(Ints)
      return Xis
end
#
function fXis_6(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(6))*H6(Ints)
      return Xis
end
#
function fXis_7(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(7))*H7(Ints)
      return Xis
end
#
function fXis_8(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(8))*H8(Ints)
      return Xis
end
#
function fXis_9(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(9))*H9(Ints)
      return Xis
end
#
function fXis_10(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(10))*H10(Ints)
      return Xis
end
#  
function fXis_11(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(11))*H11(Ints)
      return Xis
end
#
function fXis_12(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(12))*H12(Ints)
      return Xis
end
#
function fXis_13(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(13))*H13(Ints)
      return Xis
end
#
function fXis_14(nb,TiempoC,BM) 
      Ints=BM
      Xis=sqrt(factorial(14))*H14(Ints)
      return Xis
end
#
function fXis_15(nb,TiempoC,BM)
      Ints=BM
      Xis=sqrt(factorial(15))*H15(Ints)
      return Xis
      end

      ####hermite polynomial
      
function H2(x)
      h=((x .^ 2) .- 1)/sqrt(2)
      return h
end
#  
function H3(x)
      h=(x .* ((x .^ 2).-3))/sqrt(factorial(3))
      return h
end
#
function H4(x)
      h=(x.^ 4 .- 6 .* (x .^ 2).+ 3)/sqrt(factorial(4))
      return h
end
#      
function H5(x)
      h=(x .*((x .^4) .-10 .*(x .^2) .+ 15))/sqrt(factorial(5))
      return h
end
#   
function H6(x)
      h=((x .^6) .-15 .*(x .^4) .+45*(x .^2) .-15)/sqrt(factorial(6))
      return h
end
#
function H7(x)
      h=(x .*((x .^6) .-21 .*(x .^4).+105 .*(x .^2) .-105))/sqrt(factorial(7))
      return h
end
#     
function H8(x)
      h=((x .^8) .-28 .*(x .^6) .+210 .*(x .^4)-420 .*(x .^2) .+105)/sqrt(factorial(8))
      return h
end
#   
function H9(x)
      h=(x .*((x .^8) .-36 .*(x .^6) .+378 .*(x .^4) .-1260 .*(x .^2) .+945))/sqrt(factorial(9))
      return h
end
#
function H10(x)
      h=((x^10)-45*(x^8)+630*(x^6)-3150*(x^4)+4725*(x^2)-945)/sqrt(factorial(10))
      return h
end
#
function H11(x)
      h=(x*(((x^10)-55*(x^8)+990*(x^6)-6930*(x^4)+1735*(x^2)-10395)))/sqrt(factorial(11))
      return h
end
# 
function H12(x)
      h=((x^12)-66*(x^10)+1485*(x^8)-13860*(x^6)+51975*(x^4)-62370*(x^2)+10395)/sqrt(factorial(12))
      return h
end
#
function H13(x)
      h=(x*((x^12)-78*(x^10)+2145*(x^8)-25740*(x^6)+135135*(x^4)-270270*(x^2)+135135))/sqrt(factorial(13))
      return h
end
#
function H14(x)
      h=((x^14)-91*(x^12)+3003*(x^10)-45045*(x^8)+315315*(x^6)-945945*(x^4)+945945*(x^2)-135135)/sqrt(factorial(14))
      return h
end
function H15(x)
      h=(x*((x^14)-105*(x^12)+4095*(x^10)-75075*(x^8)+675675*(x^6)-2837835*(x^4)+4729725*(x^2)-2027025))/sqrt(factorial(15))
      return h
end
#
function fXis_2_j(js,TiempoC,BM)
      Ints=BM
      Xis=Ints[js[1]]*Ints[js[2]]
      return Xis
end
#
function fXis_3_12(js,TiempoC,BM)
      Ints=BM
      Xis=Ints[js[1]]*H2(Ints[js[2]])
      return Xis
end
#  
function fXis_3_21(js,TiempoC,BM)
      Ints=BM
      Xis=sqrt(2)*H2(Ints[js[1]])*Ints[js[2]]
      return Xis
end