/**
 * Created by jon on 02/10/15.
 */

$(document).ready(function() {
    $(".action").hover(
        function(event) {
            if (!$(this).hasClass("exempt")) {
                $(this).toggleClass("hovering");
            }
        }
    );

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
});