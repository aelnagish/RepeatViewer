
(function(){
  var box=document.getElementById('search'),out=document.getElementById('results');
  if(!box||!window.MANUAL_INDEX)return;
  function esc(s){return s.replace(/[&<>]/g,function(c){return{'&':'&amp;','<':'&lt;','>':'&gt;'}[c]})}
  box.addEventListener('input',function(){
    var q=box.value.trim().toLowerCase();
    if(q.length<2){out.hidden=true;out.innerHTML='';return}
    var words=q.split(/\s+/),hits=[];
    MANUAL_INDEX.forEach(function(s){
      var t=(s.t+' '+s.x).toLowerCase(),ok=words.every(function(w){return t.indexOf(w)>=0});
      if(ok){var i=s.x.toLowerCase().indexOf(words[0]);
        hits.push({s:s,snip:i<0?s.x.slice(0,140):s.x.slice(Math.max(0,i-60),i+100)})}
    });
    out.hidden=false;
    out.innerHTML='<h4>'+hits.length+' section'+(hits.length==1?'':'s')+' found</h4><ol>'+hits.slice(0,30).map(function(h){
      return '<li><a href="'+h.s.u+'">'+esc(h.s.t)+'</a><p>…'+esc(h.snip)+'…</p></li>'}).join('')+'</ol>';
  });
  document.addEventListener('keydown',function(e){if(e.key=='/'&&document.activeElement!==box){e.preventDefault();box.focus()}});
})();
