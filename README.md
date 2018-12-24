gitの使い方が良くわかってないので変なところがあっても気にしないでください

# dfttest
181224<br>
-ベースをdffttest1.9.4.3に変更<br>
-SSE/AVX/AVX2の最適化調整 あくまで自分の環境(i7-4790K Haswell)で速くなるようにしただけ (AVX2はあまり効果がないが)<br>
 opt=4をAVX2用にした<br>
-アライメント調整したらちょっと早くなったような気がする<br>
-dither処理最適化 (ditherの値↓により処理分岐させるようにした)<br>
 dither<100<br>
  1スレッド用最適化コード<br>
 100<=dither<200<br>
  dither処理が律速になるとき用にdll内でスレッド分割 <br>
  (逆に、avisynth+でMTを使うときはdither<100の方が無意味にスレッドが増えなくて良い)<br>
 200<=dither<300<br>
  従来互換 元々のコード(なので遅い)<br>
<br>
-バイナリはSSE2でコンパイルされているので、SSE2以上必須<br>
 openmpを使っているのでvcomp140.dllが要る(VisualC++2015/2017再配布用パッケージをインストールすること)<br>
 x64はちゃんと確認してない(一応コンパイルが通るようにはしてあるだけ)<br>
 <br>
<br>
this is fork of dfttest 1.9.4.3 and some speed tune for my enviroment is added.<br>
mainly tuned at dither function(it costs wastefully a lot).<br>
And not tested a lot.<br>

