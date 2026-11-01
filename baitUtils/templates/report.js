// Table of contents: highlight the section in view. Everything else in the
// report works without JavaScript (native <details> for collapsible blocks).
(function () {
  var links = Array.prototype.slice.call(document.querySelectorAll('.toc a[href^="#"]'));
  if (!links.length || !('IntersectionObserver' in window)) { return; }
  var byId = {};
  links.forEach(function (a) { byId[a.getAttribute('href').slice(1)] = a; });
  var observer = new IntersectionObserver(function (entries) {
    entries.forEach(function (entry) {
      if (entry.isIntersecting) {
        links.forEach(function (a) { a.classList.remove('active'); });
        var link = byId[entry.target.id];
        if (link) { link.classList.add('active'); }
      }
    });
  }, { rootMargin: '-10% 0px -80% 0px' });
  Object.keys(byId).forEach(function (id) {
    var el = document.getElementById(id);
    if (el) { observer.observe(el); }
  });
})();
