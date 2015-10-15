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
});