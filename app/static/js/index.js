/**
 * Created by jon on 02/10/15.
 */

$(document).ready(function() {
    $(".action").hover(
        function(event) {
            $(this).toggleClass("hovering");
        }
    );
});