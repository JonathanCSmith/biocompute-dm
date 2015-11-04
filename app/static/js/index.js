/**
 * Created by jon on 02/10/15.
 */

function addScrollListeners() {
    var didScroll;
    var lastScrollTop = 0;
    var deltaScroll = 5;
    var bigtop = $(".bigtop");
    var height = bigtop.outerHeight();

    $(window).scroll(function(event) {
        didScroll = true;
    });

    setInterval(function() {
        if (didScroll) {
            hasScrolled();
            didScroll = false;
        }
    }, 250);

    function hasScrolled() {
        var tp = $(this).scrollTop();

        if (Math.abs(lastScrollTop - tp) <= deltaScroll) {
            return;
        }

        if (tp > lastScrollTop && tp > height && bigtop.is(":visible")) {
            bigtop.hide();
            window.scrollTo(0, 0);
            lastScrollTop = 0;
            return;
        }

        else if (tp - lastScrollTop < 0 && tp < 20) {
            bigtop.show();
        }

        //else if (tp + $(window).height() < $(document).height()) {
        //    $(".bigtop").show();
        //}

        lastScrollTop = tp;
    }
}

function addStickyLocation() {
    // clean up arguments (allows prefix to be optional - a bit of overkill)
    var prefix = window.location.host.split('_').join('');
    var pre_name;

    $(".remember").click(
        function(event) {
            if (typeof(pre_name) == 'undefined') {
                pre_name = window.name;
            }

            // get the scroll positions
            var top = 0, left = 0;
            if (typeof(window.pageYOffset) == 'number') { // netscape
                top = window.pageYOffset;
                left = window.pageXOffset;
            }

            else if (document.body && (document.body.scrollLeft || document.body.scrollTop)) { // dom
                top = document.body.scrollTop;
                left = document.body.scrollLeft;
            }

            else if (document.documentElement && (document.documentElement.scrollLeft || document.documentElement.scrollTop)) { // ie6
                top = document.documentElement.scrollTop;
                left = document.documentElement.scrollLeft;
            }

            // store the scroll
            if (top || left) {
                window.name = prefix + '_' + left + '_' + top + '_' + pre_name;
            }

            return true;
        }
    );

    if (window.name.search('^'+prefix+'_(\\d+)_(\\d+)_') == 0) {
        var name = window.name.split('_');
        window.scrollTo(name[1], name[2]);
        window.name = name.slice(3).join('_');
    }
}

function addHovers() {
    $(".action").hover(
        function(event) {
            if (!$(this).hasClass("exempt")) {
                $(this).toggleClass("hovering");
            }
        }
    );
}

$(document).ready(function() {
    addScrollListeners();
    addStickyLocation();
    addHovers();
});