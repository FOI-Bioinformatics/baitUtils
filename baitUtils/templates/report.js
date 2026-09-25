        // Collapsible sections
        document.querySelectorAll('.collapsible-header').forEach(function(header) {
            header.addEventListener('click', function() {
                this.classList.toggle('active');
                var content = this.nextElementSibling;
                content.classList.toggle('active');
                
                if (content.classList.contains('active')) {
                    content.style.maxHeight = content.scrollHeight + 'px';
                } else {
                    content.style.maxHeight = '0';
                }
            });
        });
        
        // Smooth scrolling for internal links
        document.querySelectorAll('a[href^="#"]').forEach(anchor => {
            anchor.addEventListener('click', function (e) {
                e.preventDefault();
                document.querySelector(this.getAttribute('href')).scrollIntoView({
                    behavior: 'smooth'
                });
            });
        });
        
