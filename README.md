gitの使い方が良くわかってないので変なところがあっても気にしないでください

# dfttest
dfttestをプログラミングの勉強を兼ねて弄ってます

-ベースをdffttest1.9.4.3に変更<br>
-SSE/AVX/AVX2の最適化調整 あくまで自分の環境(i7-4790K Haswell)で速くなるようにしただけ (AVX2はあまり効果がないが)<br>
 opt=4をAVX2用にした<br>

this is fork of dfttest 1.9.4.3 and some speed tune for my enviroment is added.<br>
mainly tuned at dither function(it costs wastefully a lot).<br>
And not tested a lot.<br>

